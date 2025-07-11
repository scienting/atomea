from typing import TYPE_CHECKING, Iterator

from abc import ABC

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class _BaseLogicalOperator(SelectionExpression, ABC):
    """Abstract base class for logical operations (AND, OR, NOT)."""

    def __init__(
        self,
        left: SelectionExpression,
        right: SelectionExpression | None = None,
        operand: SelectionExpression | None = None,
    ):
        """
        Args:
            left: The `SelectionExpression` on the left side of the logical operation.
            right: The `SelectionExpression` on the right side of the logical operation.
            operand: The `SelectionExpression` to negate.
        """
        self._left = left
        self._right = right
        self._operand = operand


class AndExpression(_BaseLogicalOperator):
    """Represents a logical AND operation between two selection expressions.

    This expression evaluates to `True` for an atom if and only if both
    the left and right sub-expressions evaluate to `True` for that atom.
    """

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool_]]:
        """Evaluates the logical AND of the left and right sub-expressions.

        For each microstate, this method retrieves the boolean masks from
        both the left and right sub-expressions and applies a NumPy
        logical AND operation element-wise.

        Args:
            ensemble: The `Ensemble` object containing the data.
            run_id: The ID of the run within the ensemble.
            micro_id: Microstate IDs to evaluate expression for. If `None`, then it
                evaluates for every single possible microstate.

        Yields:
            A 1D boolean NumPy array (`npt.NDArray[np.bool_]`) for each microstate,
                representing the element-wise logical AND of the masks produced by
                the left and right sub-expressions.
        """
        left_masks = self._left.evaluate(ensemble, run_id, micro_id)
        right_masks = self._right.evaluate(ensemble, run_id, micro_id)
        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_and(l_mask, r_mask)


class OrExpression(_BaseLogicalOperator):
    """Represents a logical OR operation between two selection expressions.

    This expression evaluates to `True` for an atom if either the left or
    the right sub-expression (or both) evaluate to `True` for that atom.
    """

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool_]]:
        """Evaluates the logical OR of the left and right sub-expressions.

        For each microstate, this method retrieves the boolean masks from
        both the left and right sub-expressions and applies a NumPy
        logical OR operation element-wise.

        Args:
            ensemble: The `Ensemble` object containing the data.
            run_id: The ID of the run within the ensemble.
            micro_id: Microstate IDs to evaluate expression for. If `None`, then it
                evaluates for every single possible microstate.

        Yields:
            A 1D boolean NumPy array (`npt.NDArray[np.bool_]`) for each microstate,
                representing the element-wise logical OR of the masks produced by
                the left and right sub-expressions.
        """
        left_masks = self._left.evaluate(ensemble, run_id, micro_id)
        right_masks = self._right.evaluate(ensemble, run_id, micro_id)
        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_or(l_mask, r_mask)


class NotExpression(_BaseLogicalOperator):
    """Represents a logical NOT operation on a selection expression.

    This expression evaluates to `True` for an atom if its sub-expression
    evaluates to `False`, and vice-versa.
    """

    def __init__(self, expression: SelectionExpression):
        super().__init__(None, None, expression)  # Use operand for unary
        self.expression = expression

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool_]]:
        """Evaluates the logical NOT of the sub-expression.

        For each microstate, this method retrieves the boolean mask from
        its sub-expression and applies a NumPy logical NOT operation element-wise.

        Args:
            ensemble: The `Ensemble` object containing the data.
            run_id: The ID of the run within the ensemble.
            micro_id: Microstate IDs to evaluate expression for. If `None`, then it
                evaluates for every single possible microstate.

        Yields:
            A 1D boolean NumPy array (`npt.NDArray[np.bool_]`) for each microstate,
                representing the element-wise logical NOT of the mask produced by
                the sub-expression.
        """
        masks = self.expression.evaluate(ensemble, run_id, micro_id)
        for mask in masks:
            yield np.logical_not(mask)


class SelectionOperator:
    """
    Provides static methods to act as factories for SelectionExpression instances.
    This is often called from from
    """

    @staticmethod
    def AND(left: SelectionExpression, right: SelectionExpression) -> AndExpression:
        """Creates an [AndExpression][selection.expressions.ops.AndExpression]."""
        if right is None:  # Added for robustness, though type hints should catch it.
            raise ValueError("AND operator requires a right-hand side expression.")
        return AndExpression(left, right)

    @staticmethod
    def OR(left: SelectionExpression, right: SelectionExpression) -> OrExpression:
        """Creates an [OrExpression][selection.expressions.ops.OrExpression]."""
        if right is None:
            raise ValueError("OR operator requires a right-hand side expression.")
        return OrExpression(left, right)

    @staticmethod
    def NOT(expression: SelectionExpression) -> NotExpression:
        """Creates a [NotExpression][selection.expressions.ops.NotExpression]."""
        return NotExpression(expression)
