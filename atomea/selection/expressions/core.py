from typing import TYPE_CHECKING, Iterator

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class SelectionExpression(ABC):
    """Abstract base class for all selection criteria and logical operations.

    This class defines the interface for expressions that can be evaluated
    to produce a boolean mask for atoms within an `Ensemble` microstate.
    It also provides dunder methods for convenient logical operations
    (AND, OR, NOT) using Python's `&`, `|`, and `~` operators.
    """

    @abstractmethod
    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool_]]:
        """Evaluates this expression for specified microstates.

        This method generates an atom-level boolean mask for each requested microstate.
        The mask indicates which atoms meet the criteria defined by the expression.

        Args:
            ensemble: The `Ensemble` object containing the data to be evaluated.
            run_id: The ID of the specific run within the ensemble to evaluate.
                If `None`, the expression might apply across all runs or use a default.
            micro_id: Specifies the microstate IDs for which to evaluate the expression.
                If `None`, the expression is evaluated for every single possible microstate
                available for the given run. This can be an integer, a slice object,
                or an iterable of integers.

        Yields:
            A 1D boolean NumPy array (`npt.NDArray[np.bool_]`). Each element in the
                array corresponds to an atom in the specified microstate, with `True`
                indicating that the atom meets the selection criteria and `False` otherwise.
                Each yielded array corresponds to one microstate.
        """

    def __and__(self, other: "SelectionExpression") -> "SelectionExpression":
        """Performs a logical AND operation with another SelectionExpression.

        Args:
            other: The other `SelectionExpression` to combine with a logical AND.

        Returns:
            An `AndExpression` instance representing the logical AND of this
            expression and the `other` expression.
        """
        return AndExpression(self, other)

    def __or__(self, other: "SelectionExpression") -> "SelectionExpression":
        """Performs a logical OR operation with another SelectionExpression.

        Args:
            other: The other `SelectionExpression` to combine with a logical OR.

        Returns:
            An `OrExpression` instance representing the logical OR of this
            expression and the `other` expression.
        """
        return OrExpression(self, other)

    def __invert__(self) -> "SelectionExpression":
        """Performs a logical NOT operation on this SelectionExpression.

        Returns:
            A `NotExpression` instance representing the logical NOT of this expression.
        """
        return NotExpression(self)


class BinaryLogicalExpression(SelectionExpression, ABC):
    """Abstract base class for binary logical operations (AND, OR).

    This class serves as a common base for expressions that combine two
    `SelectionExpression` instances using a binary logical operator.
    """

    def __init__(self, left: SelectionExpression, right: SelectionExpression):
        """Initializes a BinaryLogicalExpression.

        Args:
            left: The `SelectionExpression` on the left side of the binary operation.
            right: The `SelectionExpression` on the right side of the binary operation.
        """
        self.left = left
        self.right = right


class AndExpression(BinaryLogicalExpression):
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
        left_masks = self.left.evaluate(ensemble, run_id, micro_id)
        right_masks = self.right.evaluate(ensemble, run_id, micro_id)

        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_and(l_mask, r_mask)


class OrExpression(BinaryLogicalExpression):
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
        left_masks = self.left.evaluate(ensemble, run_id, micro_id)
        right_masks = self.right.evaluate(ensemble, run_id, micro_id)

        for l_mask, r_mask in zip(left_masks, right_masks):
            yield np.logical_or(l_mask, r_mask)


class NotExpression(SelectionExpression):
    """Represents a logical NOT operation on a selection expression.

    This expression evaluates to `True` for an atom if its sub-expression
    evaluates to `False`, and vice-versa.
    """

    def __init__(self, expression: SelectionExpression):
        """Initializes a NotExpression.

        Args:
            expression: The `SelectionExpression` to apply the logical NOT operation to.
        """
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
