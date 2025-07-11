from typing import TYPE_CHECKING, Iterator

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from loguru import logger

import atomea.typing as adt
from atomea.data import OptionalSliceSpec

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class SelectionExpression(ABC):
    """Abstract base class for all selection criteria and logical operations.

    This class defines the interface for expressions that can be evaluated
    to produce a boolean mask for atoms within an `Ensemble` microstate.
    It also provides dunder methods for convenient logical operations
    (`AND`, `OR`, `NOT`) using Python's `&`, `|`, and `~` operators.
    """

    def init_mask(self) -> adt.Bool:
        logger.warning("Did not find any atom types; eval will be false.")
        num_atoms: int | None = self._n_atoms if self._n_atoms is not None else 0
        mask: adt.Bool = np.zeros(num_atoms, dtype=np.dtype(np.bool))
        return mask

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
