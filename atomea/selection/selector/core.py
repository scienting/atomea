from typing import TYPE_CHECKING

import numpy as np

from atomea.selection.expressions import (
    AndExpression,
    AtomTypeIs,
    DistanceWithin,
    MolIdIs,
    NotExpression,
    OrExpression,
    SelectionExpression,
)

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class EnsembleSelector:
    """
    A fluent API for building complex selection queries on an Ensemble's data.

    Allows chaining of selection criteria with logical AND, OR, and NOT operations.
    The selection is lazy, meaning the actual computation happens only when
    `get_mask()` is called.
    """

    def __init__(self, ensemble: "Ensemble", run_id: str | None = None):
        """
        Initializes the EnsembleSelector for a specific ensemble and run.

        Args:
            ensemble: The Ensemble object to perform selections on.
            run_id: The ID of the run within the ensemble for which to generate masks.
        """
        self._ensemble = ensemble
        self._run_id = run_id
        self._current_expression: SelectionExpression | None = None
        # This flag helps manage implicit AND vs explicit OR/AND
        self._next_op_is_or: bool = False

    def _add_expression(self, new_expr: SelectionExpression) -> "EnsembleSelector":
        """Internal method to add a new expression to the tree."""
        if self._current_expression is None:
            self._current_expression = new_expr
        else:
            if self._next_op_is_or:
                self._current_expression = OrExpression(
                    self._current_expression, new_expr
                )
                self._next_op_is_or = False  # Reset after applying OR
            else:
                # Default is AND
                self._current_expression = AndExpression(
                    self._current_expression, new_expr
                )
        return self

    def mol_id_is(self, mol_ids: list[int]) -> "EnsembleSelector":
        """
        Selects atoms whose molecule ID is present in the provided list.

        Args:
            mol_ids: A list of integer molecule IDs to select.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        return self._add_expression(MolIdIs(mol_ids))

    def atom_type_is(self, atom_types: list[str]) -> "EnsembleSelector":
        """
        Selects atoms whose atom type is present in the provided list.

        Args:
            atom_types: A list of string atom types (e.g., "C", "N", "H") to select.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        return self._add_expression(AtomTypeIs(atom_types))

    def distance_within(
        self, center_atom_index: int, radius: float
    ) -> "EnsembleSelector":
        """
        Selects atoms within a specified radial distance from a given center atom index.
        This filter is applied per microstate.

        Args:
            center_atom_index: The 0-based index of the atom to use as the center.
            radius: The maximum distance (inclusive) from the center atom in the same units as coordinates.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        return self._add_expression(DistanceWithin(center_atom_index, radius))

    def and_(self) -> "EnsembleSelector":
        """
        Sets the logical operator for the next chained filter to AND.
        This is the default behavior if no logical operator is explicitly called.
        """
        # No explicit action needed here as _add_expression defaults to AND.
        # This method primarily serves to make the chaining explicit and readable.
        self._next_op_is_or = False  # Ensure it's explicitly AND
        return self

    def or_(self) -> "EnsembleSelector":
        """
        Sets the logical operator for the next chained filter to OR.
        The next filter added will be OR'd with the current expression tree.
        """
        if self._current_expression is None:
            raise ValueError(
                "Cannot start a query with an 'or_' operator. Add a filter first."
            )
        self._next_op_is_or = True
        return self

    def not_(self) -> "EnsembleSelector":
        """
        Applies a logical NOT operation to the *entire* current expression tree.
        This should typically be used at the end of a sub-expression or for a single filter.
        Example: `select().not_().mol_id_in([0])`
        Example: `select().mol_id_in([0]).and_().not_().atom_type_is(["C"])`
        """
        if self._current_expression is None:
            raise ValueError(
                "Cannot apply 'not_' to an empty selection. Add a filter first."
            )
        self._current_expression = NotExpression(self._current_expression)
        return self

    def get_mask(self) -> dict[int, np.ndarray]:
        """
        Executes the built selection query and returns a dictionary of atom-level
        boolean masks, one for each microstate.

        The mask for each microstate will have a length equal to the number of atoms
        in that microstate, indicating which atoms meet the criteria.

        Returns:
            A dictionary where keys are `microstate_id` (int) and values are
            1D boolean NumPy arrays representing the atom selection mask for that microstate.
        """
        if self._current_expression is None:
            # If no filters were added, return all True for all atoms in all microstates
            all_microstate_ids = self._get_all_microstate_ids()
            if not all_microstate_ids:
                return {}

            masks = {}
            for ms_id in all_microstate_ids:
                coords = self._ensemble.coordinates.read(
                    ens_id=self._ensemble.label,
                    run_id=self._run_id,
                    microstate_id=ms_id,
                )
                num_atoms = coords.shape[1] if coords is not None else 0
                masks[ms_id] = np.ones(num_atoms, dtype=bool)
            return masks

        # Get all microstate IDs that exist for this ensemble and run
        all_microstate_ids = self._get_all_microstate_ids()
        if not all_microstate_ids:
            return {}

        result_masks: dict[int, np.ndarray] = {}
        for microstate_id in all_microstate_ids:
            # Evaluate the entire expression tree for each microstate
            mask = self._current_expression.evaluate(
                self._ensemble, microstate_id, self._run_id
            )
            result_masks[microstate_id] = mask
        return result_masks

    def _get_all_microstate_ids(self) -> list[int]:
        """Helper to get all unique microstate IDs for the current ensemble and run."""
        # Assuming prj.energy.potential_mm is a reliable source for all microstate IDs
        energy_df = self._ensemble._parent.energy.potential_mm.read(
            ens_id=self._ensemble.label, run_id=self._run_id
        )
        if energy_df is None or energy_df.shape[0] == 0:
            # If no energy data, try to get microstate IDs from coordinates
            # This is less reliable as coordinates might not exist for all microstates
            # or might be sparse.
            # A more robust solution might involve a dedicated 'microstates' table.
            coords_data = self._ensemble.coordinates.get(
                ens_id=self._ensemble.label, run_id=self._run_id
            )
            if coords_data is not None:
                # If Zarr array, we might need to inspect its structure for microstate dimension
                # This is a simplification; actual implementation might need to query Zarr metadata
                # For now, assume it's just the microstate index.
                # This part is highly dependent on how Zarr is structured for microstates.
                # If Zarr stores (microstate_idx, num_atoms, 3), then coords_data.shape[0]
                # would give the number of microstates.
                # For now, let's assume we can get a list of available microstate IDs from the store.
                # This needs to be refined based on actual Zarr/Polars usage.
                # Placeholder:
                # If coords_data is a Zarr array, its shape[0] would be the number of microstates
                # Assuming microstate IDs are 0 to N-1
                return list(range(coords_data.shape[0]))
            return []  # No microstates found

        return sorted(energy_df["microstate_id"].unique().to_list())
