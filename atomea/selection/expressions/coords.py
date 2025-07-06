from typing import TYPE_CHECKING, Iterator

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec, SliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class DistanceWithin(SelectionExpression):
    """
    Selects atoms within a certain distance from a reference selection.
    Selections include the reference atoms.
    """

    def __init__(self, from_atoms: SliceSpec, dist: float):
        """
        Args:
            atoms: Reference atoms to compute distances from.
            dist: Selection distance
        """
        self.from_atoms = from_atoms
        self.dist = dist

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        if micro_id:
            view = (micro_id, slice(None), slice(None))
        else:
            view = (slice(None),)

        for chunk in ensemble.coordinates.iter(run_id=run_id, view=view, chunk_size=1):
            if chunk is None or chunk.shape[1] == 0:
                # No coordinates for this microstate, return empty mask
                yield np.array([], dtype=np.dtype(np.bool))
            else:
                from_coords = chunk[self.from_atoms]
                distances = np.linalg.norm(chunk - from_coords, axis=1)
                yield distances <= self.dist
