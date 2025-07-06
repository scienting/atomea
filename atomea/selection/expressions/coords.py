from typing import TYPE_CHECKING

import numpy as np

from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class DistanceWithin(SelectionExpression):
    """
    Selects atoms within a certain distance from a center atom.
    This is a MICROSTATE cadence filter.
    """

    def __init__(self, center_atom_index: int, radius: float):
        self.center_atom_index = center_atom_index
        self.radius = radius

    def evaluate(
        self, ensemble: "Ensemble", run_id: str
    ) -> np.ndarray:
        # Coordinates are MICROSTATE cadence, so load per microstate
        coords = ensemble.coordinates.values(
            run_id=run_id
        )
        if coords is None or coords.shape[1] == 0:
            # No coordinates for this microstate, return empty mask
            return np.array([], dtype=bool)

        # Assuming coords is (1, num_atoms, 3) or (num_atoms, 3) for a single microstate
        # We need to handle the case where `values` returns (num_microstates, num_atoms, 3)
        # For evaluate, we expect data for *one* microstate.
        if coords.ndim == 3 and coords.shape[0] == 1:
            coords_single_microstate = coords[0]
        elif coords.ndim == 2:
            coords_single_microstate = coords
        else:
            raise ValueError(
                f"Unexpected coordinate shape for single microstate: {coords.shape}"
            )

        if self.center_atom_index >= coords_single_microstate.shape[0]:
            raise IndexError(
                f"Center atom index {self.center_atom_index} out of bounds for {coords_single_microstate.shape[0]} atoms."
            )

        center_coords = coords_single_microstate[self.center_atom_index]
        distances = np.linalg.norm(coords_single_microstate - center_coords, axis=1)
        return distances <= self.radius
