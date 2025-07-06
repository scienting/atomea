from typing import TYPE_CHECKING, Iterator

import numpy as np
import numpy.typing as npt

from atomea.data import SliceSpec, OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class DistanceWithin(SelectionExpression):
    """
    Selects atoms within a certain distance from a selection.
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
            view = (slice(None), slice(None))
        else:
            view = 

        for coords in ensemble.coordinates.iter(run_id, 
        coords = ensemble.coordinates.read(view=view, run_id=run_id)
        if coords is None or coords.shape[1] == 0:
            # No coordinates for this microstate, return empty mask
            return np.array([], dtype=np.dtype(np.bool))

        for coord in coords:
            from_coords = coord[self.from_atoms]
            distances = np.linalg.norm(coord - from_coords, axis=1)
            return distances <= self.dist
