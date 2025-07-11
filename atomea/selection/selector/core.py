from typing import TYPE_CHECKING, Callable, Iterator, Self

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec, SliceSpec
from atomea.selection.expressions import (
    AtomTypeIs,
    DistanceWithin,
    MolIdIs,
    SelectionExpression,
    SelectionOperator,
)

# See https://stackoverflow.com/questions/33533148/how-do-i-type-hint-a-method-with-the-type-of-the-enclosing-class


if TYPE_CHECKING:
    from atomea.containers import Ensemble


class EnsembleSelector:
    """
    A fluent API for building complex selection queries on an Ensemble's data.

    Allows chaining of selection criteria with logical `.AND()`, `.OR()`, and `.NOT()` operations.
    The selection is lazy, meaning the actual computation happens only when
    `get_mask()` is called.
    """

    def __init__(
        self, ensemble: "Ensemble", run_id: str | None = None, is_negated: bool = False
    ):
        """
        Initializes the EnsembleSelector for a specific ensemble and run.

        Args:
            ensemble: The Ensemble object to perform selections on.
            run_id: The ID of the run within the ensemble for which to generate masks.
            is_negated: Any filter after this specific instance will be negated.
        """
        self._ensemble = ensemble
        """Ensemble object to perform selections on."""

        self._run_id = run_id
        """Run ID to use if needed for expressions."""

        self._current_expression: SelectionExpression | None = None
        """This flag helps manage implicit AND vs explicit OR/AND."""

        self._pending_operator_factory: Callable[
            [SelectionExpression, SelectionExpression | None], SelectionExpression
        ] = SelectionOperator.AND
        """Indicates how to apply the next operation."""

        self._is_negated = is_negated
        """Indicates if we should negate the next selection by using `NOT`"""

        self._n_atoms: int | None = ensemble.topology.atoms.n_atoms
        """Number of atoms to initialize the boolean mask with."""

    def _add_expression(self, new_expr: SelectionExpression) -> "EnsembleSelector":
        """Internal method to add a new expression to the tree."""
        # Apply negation immediately if this selector instance is marked as negated
        if self._is_negated:
            new_expr = SelectionOperator.NOT(new_expr)
            # IMPORTANT: After applying the negation, create a *new* non-negated selector
            # for subsequent chaining. This prevents `NOT().NOT().atom_types()` from negating twice.
            new_selector = EnsembleSelector(self._ensemble, self._run_id)
            new_selector._current_expression = new_expr
            new_selector._pending_operator_factory = self._pending_operator_factory
            return new_selector

        if self._current_expression is None:
            self._current_expression = new_expr
        else:
            self._current_expression = self._pending_operator_factory(
                self._current_expression, new_expr
            )

        if self._pending_operator_factory is not SelectionOperator.OR:
            self._pending_operator_factory = SelectionOperator.AND

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
            if from_atoms._current_expression is None:
                raise ValueError(
                    "The 'from_atoms' selector has no active selection. Add filters to it first."
                )
            # Ensure the nested selector's expression is processed correctly
            nested_expr = from_atoms._current_expression
            return self._add_expression(DistanceWithin(nested_expr, dist))
        else:
            return self._add_expression(DistanceWithin(from_atoms, dist))

    def AND(self) -> "EnsembleSelector":
        """
        Sets the logical operator for the next chained filter to AND.
        This is the default behavior if no logical operator is explicitly called.
        """
        self._pending_operator_factory = SelectionOperator.AND
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
        self._pending_operator_factory = SelectionOperator.OR
        return self

    def NOT(self) -> "EnsembleSelector":
        """
        Returns a _new_ EnsembleSelector instance where the next added filter will be
        negated. This handles `select().NOT().atom_types(...)`.
        """
        # Create a new selector instance that is marked as negated.
        # This allows chaining like `select().NOT().atom_types(...)`
        # or `select().atom_types(...).AND().NOT().molecule_ids(...)`
        negated_selector = EnsembleSelector(
            self._ensemble, self._run_id, is_negated=True
        )
        # Carry over the current expression and pending operator from the original selector
        negated_selector._current_expression = self._current_expression
        negated_selector._pending_operator_factory = self._pending_operator_factory
        return negated_selector

    def get_mask(
        self, micro_id: OptionalSliceSpec = None
    ) -> Iterator[npt.NDArray[np.bool]]:
        """
        Executes the built selection query and yields atom-level
        boolean masks, one for each microstate.

        The mask for each microstate will have a length equal to the number of atoms in that microstate, indicating which atoms meet the criteria.

        Args:
            micro_id: Specifies the microstate IDs for which to generate masks.
                Can be an integer, a slice, or None (to get masks for all microstates).

        Returns:
            1D boolean NumPy arrays representing the atom selection mask for that microstate.
        """
        if self._current_expression is None:
            # If no selection expression, yield all-true masks for the specified microstates.
            for _ in range(
                0, self._ensemble.n_micro(view=micro_id, run_id=self._run_id)
            ):
                num_atoms = self._n_atoms if self._n_atoms is not None else 0
                yield np.ones(num_atoms, dtype=bool)
        else:
            # If there is a selection expression, delegate the evaluation directly to it.
            # The SelectionExpression.evaluate method already returns an iterator of masks.
            yield from self._current_expression.evaluate(
                self._ensemble, self._run_id, micro_id=micro_id
            )
