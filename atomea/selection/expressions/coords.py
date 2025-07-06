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

    def __init__(self, from_atoms: SelectionExpression | SliceSpec, dist: float):
        """
        Args:
            atoms: Reference atoms to compute distances from. This can be
                a `SliceSpec` or another `SelectionExpression`, in which
                case the distance will be computed from the atoms selected
                by that expression.
            dist: Selection distance
        """
        self.from_atoms = from_atoms
        self.dist = dist

    def _evaluate_selection(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        reference_masks_iterator = self.from_atoms.evaluate(ensemble, run_id, micro_id)

        # We need to iterate through both the coordinates and the reference masks
        # to ensure they are synchronized per microstate.
        for (ms_coords, ms_id), ref_mask in zip(
            ensemble.coordinates.iter(
                run_id=run_id,
                view=(micro_id, slice(None), slice(None)),
                return_microstate_id=True,
            ),
            reference_masks_iterator,
        ):
            if ms_coords is None or ms_coords.shape[1] == 0:
                yield np.array([], dtype=np.dtype(np.bool))
                continue

            if not np.any(ref_mask):  # If no atoms are selected by the reference mask
                yield np.zeros(ms_coords.shape[0], dtype=np.dtype(np.bool))
                continue

            # Get coordinates of reference atoms based on the mask
            from_coords = ms_coords[ref_mask]

            # Compute distances from all atoms to the closest reference atom
            # Reshape for broadcasting: (num_all_atoms, 1, 3) - (1, num_ref_atoms, 3)
            # This calculates all pairwise distances
            diff = ms_coords[:, np.newaxis, :] - from_coords[np.newaxis, :, :]
            distances = np.linalg.norm(diff, axis=2)  # (num_all_atoms, num_ref_atoms)

            # Find the minimum distance from any reference atom for each atom
            min_distances = np.min(distances, axis=1)

            yield min_distances <= self.dist

    def _evaluate_indices(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        # Original logic for direct atom indices/slices
        view = (micro_id, slice(None), slice(None)) if micro_id else (slice(None),)

        for chunk in ensemble.coordinates.iter(run_id=run_id, view=view, chunk_size=1):
            if chunk is None or chunk.shape[1] == 0:
                yield np.array([], dtype=np.dtype(np.bool))
                continue

            # Ensure from_atoms is handled correctly for indexing
            # If from_atoms is an int, make it a list for consistent indexing
            if isinstance(self.from_atoms, int):
                from_coords = chunk[self.from_atoms : self.from_atoms + 1]
            else:  # Assume slice or list[int]
                from_coords = chunk[self.from_atoms]

            # Ensure from_coords is not empty (e.g., if a slice results in no atoms)
            if from_coords.shape[0] == 0:
                yield np.zeros(chunk.shape[0], dtype=np.dtype(np.bool))
                continue

            # Compute distances from all atoms to the closest atom in from_coords
            diff = chunk[:, np.newaxis, :] - from_coords[np.newaxis, :, :]
            distances = np.linalg.norm(diff, axis=2)
            min_distances = np.min(distances, axis=1)

            yield min_distances <= self.dist

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        if isinstance(self.from_atoms, SelectionExpression):
            return self._evaluate_selection(ensemble, run_id, micro_id)
        else:
            return self._evaluate_indices(ensemble, run_id, micro_id)
