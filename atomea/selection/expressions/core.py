from typing import TYPE_CHECKING

from abc import ABC, abstractmethod

import numpy as np

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class SelectionExpression(ABC):
    """
    Abstract base class for all selection criteria and logical operations.
    Each expression can be evaluated to produce a boolean mask for a given microstate.
    """

    @abstractmethod
    def evaluate(
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        """
        Evaluates this expression for a specific microstate and returns an atom-level boolean mask.

        Args:
            ensemble: The Ensemble object containing the data.
            run_id: The ID of the run within the ensemble.
            microstate_id: The ID of the microstate for which to generate the mask.

        Returns:
            A 1D boolean NumPy array, where each element corresponds to an atom
            in the specified microstate, indicating if it meets the selection criteria.
        """

    def __and__(self, other: "SelectionExpression") -> "SelectionExpression":
        return AndExpression(self, other)

    def __or__(self, other: "SelectionExpression") -> "SelectionExpression":
        return OrExpression(self, other)

    def __invert__(self) -> "SelectionExpression":
        return NotExpression(self)


class BinaryLogicalExpression(SelectionExpression, ABC):
    """Abstract base class for binary logical operations (AND, OR)."""

    def __init__(self, left: SelectionExpression, right: SelectionExpression):
        self.left = left
        self.right = right


class AndExpression(BinaryLogicalExpression):
    """Represents a logical AND operation between two selection expressions."""

    def evaluate(
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        left_mask = self.left.evaluate(ensemble, run_id)
        right_mask = self.right.evaluate(ensemble, run_id)
        return np.logical_and(left_mask, right_mask)


class OrExpression(BinaryLogicalExpression):
    """Represents a logical OR operation between two selection expressions."""

    def evaluate(
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        left_mask = self.left.evaluate(ensemble, run_id)
        right_mask = self.right.evaluate(ensemble, run_id)
        return np.logical_or(left_mask, right_mask)


class NotExpression(SelectionExpression):
    """Represents a logical NOT operation on a selection expression."""

    def __init__(self, expression: SelectionExpression):
        self.expression = expression

    def evaluate(
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        mask = self.expression.evaluate(ensemble, run_id)
        return np.logical_not(mask)
