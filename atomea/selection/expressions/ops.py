from typing import TYPE_CHECKING, Iterator

from abc import ABC, abstractmethod

import numpy as np

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class _BaseSelectionOperator(SelectionExpression, ABC):
    """Abstract base class for all selection operators."""

    @abstractmethod
    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
        """Evaluates a logical operator.

        Args:
            ensemble: The `Ensemble` object containing the data.
            run_id: The ID of the run within the ensemble.
            micro_id: Microstate IDs to evaluate expression for. If `None`, then it
                evaluates for every single possible microstate.

        Yields:
            A 1D boolean NumPy array for each microstate.
        """


class _BinaryLogicalOperator(_BaseSelectionOperator, ABC):
    """Abstract base class for binary logical operations (AND, OR)."""

    _left: SelectionExpression
    _right: SelectionExpression

    def __init__(self, left: SelectionExpression, right: SelectionExpression) -> None:
        """
        Args:
            left: The `SelectionExpression` on the left side of the logical operation.
            right: The `SelectionExpression` on the right side of the logical operation.
        """
        self._left = left
        self._right = right


class AndExpression(_BinaryLogicalOperator):
    """Represents a logical AND operation between two selection expressions.

    This expression evaluates to `True` for an atom if and only if both
    the left and right sub-expressions evaluate to `True` for that atom.
    """

    _left: SelectionExpression
    _right: SelectionExpression

    def __init__(
        self,
        left: SelectionExpression,
        right: SelectionExpression,
    ) -> None:
        """
        Args:
            left: The `SelectionExpression` on the left side of the logical operation.
            right: The `SelectionExpression` on the right side of the logical operation.
        """
        super().__init__(left, right)

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
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
            A 1D boolean NumPy array for each microstate, representing the
                element-wise logical AND of the masks produced by
                the left and right sub-expressions.
        """
        left_masks: Iterator[adt.Bool] = self._left.evaluate(ensemble, run_id, micro_id)
        right_masks: Iterator[adt.Bool] = self._right.evaluate(
            ensemble, run_id, micro_id
        )
        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_and(l_mask, r_mask)


class OrExpression(_BinaryLogicalOperator):
    """Represents a logical OR operation between two selection expressions.

    This expression evaluates to `True` for an atom if either the left or
    the right sub-expression (or both) evaluate to `True` for that atom.
    """

    _left: SelectionExpression
    _right: SelectionExpression

    def __init__(
        self,
        left: SelectionExpression,
        right: SelectionExpression,
    ) -> None:
        """
        Args:
            left: The `SelectionExpression` on the left side of the logical operation.
            right: The `SelectionExpression` on the right side of the logical operation.
            operand: The `SelectionExpression` to negate.
        """
        super().__init__(left, right)

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
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
        left_masks: Iterator[adt.Bool] = self._left.evaluate(ensemble, run_id, micro_id)
        right_masks: Iterator[adt.Bool] = self._right.evaluate(
            ensemble, run_id, micro_id
        )
        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_or(l_mask, r_mask)


class _UnaryLogicalOperator(_BaseSelectionOperator, ABC):
    """Abstract base class for unary logical operations (NOT)."""

    _operand: SelectionExpression

    def __init__(self, operand: SelectionExpression) -> None:
        """
        Args:
            operand: The `SelectionExpression` to negate.
        """
        self._operand = operand


class NotExpression(_UnaryLogicalOperator):
    """Represents a logical NOT operation on a selection expression.

    This expression evaluates to `True` for an atom if its sub-expression
    evaluates to `False`, and vice-versa.
    """

    _operand: SelectionExpression

    def __init__(self, operand: SelectionExpression) -> None:
        super().__init__(operand)

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
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
        masks: Iterator[adt.Bool] = self._operand.evaluate(ensemble, run_id, micro_id)
        for mask in masks:
            yield np.logical_not(mask)


class SelectionOperator:
    """
    Provides static methods to act as factories for SelectionExpression instances.
    """

    @staticmethod
    def AND(left: SelectionExpression, right: SelectionExpression) -> AndExpression:
        """Creates an [AndExpression][selection.expressions.ops.AndExpression]."""
        return AndExpression(left, right)

    @staticmethod
    def OR(left: SelectionExpression, right: SelectionExpression) -> OrExpression:
        """Creates an [OrExpression][selection.expressions.ops.OrExpression]."""
        return OrExpression(left, right)

    @staticmethod
    def NOT(operand: SelectionExpression) -> NotExpression:
        """Creates a [NotExpression][selection.expressions.ops.NotExpression]."""
        return NotExpression(operand)
