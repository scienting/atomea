from typing import TYPE_CHECKING, Iterator, Self

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec, SliceSpec
from atomea.selection.expressions import (
    AndExpression,
    AtomTypeIs,
    DistanceWithin,
    MolIdIs,
    NotExpression,
    OrExpression,
    SelectionExpression,
)

# See https://stackoverflow.com/questions/33533148/how-do-i-type-hint-a-method-with-the-type-of-the-enclosing-class


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
        """Ensemble object to perform selections on."""

        self._run_id = run_id
        """Run ID to use if needed for expressions."""

        self._current_expression: SelectionExpression | None = None
        """This flag helps manage implicit AND vs explicit OR/AND."""

        self._next_op_is_or: bool = False
        """Indicates if the next expression to be added should be negated."""

        self._next_expression_is_negated: bool = False
        """Indicates if we should negate the next selection by using `NOT`"""

    def _add_expression(self, new_expr: SelectionExpression) -> "EnsembleSelector":
        """Internal method to add a new expression to the tree."""
        if self._next_expression_is_negated:
            new_expr = NotExpression(new_expr)
            self._next_expression_is_negated = False  # Reset the flag after use

        if self._current_expression is None:
            self._current_expression = new_expr
        else:
            if self._next_op_is_or:
                self._current_expression = OrExpression(
                    self._current_expression, new_expr
                )
                self._next_op_is_or = False  # Reset after applying OR
            else:
                # Default is AND if no OR was specified
                self._current_expression = AndExpression(
                    self._current_expression, new_expr
                )
        return self

    def molecule_ids(self, mol_ids: list[int]) -> "EnsembleSelector":
        """
        Selects atoms whose molecule ID is present in the provided list.

        Args:
            mol_ids: A list of integer molecule IDs to select.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        return self._add_expression(MolIdIs(mol_ids))

    def atom_types(self, atom_types: list[str]) -> "EnsembleSelector":
        """
        Selects atoms whose atom type is present in the provided list.

        Args:
            atom_types: A list of string atom types (e.g., "C", "N", "H") to select.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        return self._add_expression(AtomTypeIs(atom_types))

    def distance_within(
        self, from_atoms: SliceSpec | Self, dist: float
    ) -> "EnsembleSelector":
        """
        Selects atoms within a specified radial distance from a given center atom index.
        This filter is applied per microstate.

        Args:
            from_atoms: The 0-based index of the atom to use as the center.
            radius: The maximum distance (inclusive) from the center atom in the same units as coordinates.

        Returns:
            The EnsembleSelector instance for chaining.
        """
        if isinstance(from_atoms, EnsembleSelector):
            # If from_atoms is an EnsembleSelector, use its current expression
            if from_atoms._current_expression is None:
                raise ValueError(
                    "The 'from_atoms' selector has no active selection. Add filters to it first."
                )
            return self._add_expression(
                DistanceWithin(from_atoms._current_expression, dist)
            )
        else:
            # Otherwise, assume it's a direct index (int)
            # The DistanceWithin class will handle SliceSpec if from_atoms is more complex
            return self._add_expression(DistanceWithin(from_atoms, dist))

    def AND(self) -> "EnsembleSelector":
        """
        Sets the logical operator for the next chained filter to AND.
        This is the default behavior if no logical operator is explicitly called.
        """
        # No explicit action needed here as _add_expression defaults to AND.
        # This method primarily serves to make the chaining explicit and readable.
        self._next_op_is_or = False
        self._next_expression_is_negated = False
        return self

    def OR(self) -> "EnsembleSelector":
        """
        Sets the logical operator for the next chained filter to OR.
        The next filter added will be OR'd with the current expression tree.
        """
        if self._current_expression is None:
            raise ValueError(
                "Cannot start a query with an 'OR' operator. Add a filter first."
            )
        self._next_op_is_or = True
        self._next_expression_is_negated = False
        return self

    def NOT(self) -> "EnsembleSelector":
        """
        Applies a logical NOT operation to the *entire* current expression tree.
        This should typically be used at the end of a sub-expression or for a single filter.
        Example: `select().NOT().mol_id_in([0])`
        Example: `select().mol_id_in([0]).AND().NOT().atom_types(["C"])`
        """
        self._next_expression_is_negated = True
        return self

    def get_mask(
        self, micro_id: OptionalSliceSpec = None
    ) -> Iterator[npt.NDArray[np.bool]]:
        """
        Executes the built selection query and yields atom-level
        boolean masks, one for each microstate.

        The mask for each microstate will have a length equal to the number of atoms
        in that microstate, indicating which atoms meet the criteria.

        Args:
            micro_id: Specifies the microstate IDs for which to generate masks.
                Can be an integer, a slice, or None (to get masks for all microstates).

        Returns:
            1D boolean NumPy arrays representing the atom selection mask for that microstate.
        """
        if self._current_expression is None:
            # If no selection expression, yield all-true masks for the specified microstates.
            # We still need to iterate through coordinates to get num_atoms for each microstate.
            microstate_data_iterator = self._ensemble.coordinates.iter(
                run_id=self._run_id,
                elements=micro_id,
                chunk_size=1,
            )
            for coords_chunk in microstate_data_iterator:
                num_atoms = coords_chunk.shape[1] if coords_chunk is not None else 0
                yield np.ones(num_atoms, dtype=bool)
        else:
            # If there is a selection expression, delegate the evaluation directly to it.
            # The SelectionExpression.evaluate method already returns an iterator of masks.
            yield from self._current_expression.evaluate(
                self._ensemble, self._run_id, micro_id=micro_id
            )
