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
        chunk_size: int = 1,
    ) -> Iterator[npt.NDArray[np.bool]]:
        reference_masks_iterator = self.from_atoms.evaluate(ensemble, run_id, micro_id)

        # We need to iterate through both the coordinates and the reference masks
        # to ensure they are synchronized per microstate.
        # 'chunk_ms_coords' will be (1, N, 3) because chunk_size=1 in the iter call below
        for chunk_ms_coords, ref_mask in zip(
            ensemble.coordinates.iter(
                run_id=run_id,
                elements=micro_id,
                view=(slice(None), slice(None)),
                return_microstate_id=True,
                chunk_size=chunk_size,  # Pass chunk_size to the iterator
            ),
            reference_masks_iterator,
        ):
            # Squeeze the first dimension (microstate dimension),
            # This transforms (1, N, 3) into (N, 3)
            ms_coords = chunk_ms_coords.squeeze(axis=0)

            # Corrected check: ms_coords.shape[0] refers to the number of atoms (N)
            if (
                ms_coords is None or ms_coords.shape[0] == 0
            ):  # Check if there are no atoms (N=0)
                yield np.array([], dtype=np.dtype(np.bool))
                continue

            if not np.any(
                ref_mask
            ):  # If no atoms are selected by the reference mask (K=0)
                # Yield a mask of all False for all N atoms in the current microstate
                yield np.zeros(ms_coords.shape[0], dtype=np.dtype(np.bool))
                continue

            # Get coordinates of reference atoms based on the mask
            # Now 'ms_coords' is (N, 3) and 'ref_mask' is (N,), so 'from_coords' will be (K, 3)
            from_coords = ms_coords[ref_mask]

            # Compute distances from all atoms to the closest reference atom
            # Reshape for broadcasting: (num_all_atoms, 1, 3) - (1, num_ref_atoms, 3)
            # This calculates all pairwise distances
            # 'ms_coords' is (N, 3), so ms_coords[:, np.newaxis, :] is (N, 1, 3)
            # 'from_coords' is (K, 3), so from_coords[np.newaxis, :, :] is (1, K, 3)
            # The difference 'diff' will correctly be (N, K, 3)
            diff = ms_coords[:, np.newaxis, :] - from_coords[np.newaxis, :, :]
            distances = np.linalg.norm(diff, axis=2)  # (N, K)

            # Find the minimum distance from any reference atom for each atom
            min_distances = np.min(distances, axis=1)  # (N,)

            yield min_distances <= self.dist

    def _evaluate_indices(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
        chunk_size: int = 1,
    ) -> Iterator[npt.NDArray[np.bool]]:
        view = (slice(None), slice(None))

        for chunk in ensemble.coordinates.iter(
            run_id=run_id, elements=micro_id, view=view, chunk_size=chunk_size
        ):
            if chunk is None or chunk.shape[1] == 0:
                yield np.array([], dtype=np.dtype(np.bool))
                continue

            # Ensure from_atoms is handled correctly for indexing
            # If from_atoms is an int, make it a list for consistent indexing
            if isinstance(self.from_atoms, int):
                from_coords = chunk[
                    :, self.from_atoms : self.from_atoms + 1, :
                ]  # from_coords is (M, 1, 3)
            else:  # Assume slice or list[int]
                from_coords = chunk[:, self.from_atoms, :]  # from_coords is (M, K, 3)

            # If 'from_atoms' selection results in no atoms (K=0),
            # then from_coords will have shape (M, 0, 3).
            # In this case, no atom can be within distance, so yield all False masks.
            if from_coords.shape[1] == 0:
                num_atoms_in_ms = chunk.shape[1]
                for _ in range(chunk.shape[0]):
                    yield np.zeros(num_atoms_in_ms, dtype=np.bool)
                continue

            # Iterate through each microstate within the current 'chunk'
            # This allows 'diff' to be 3D for each individual microstate
            for i_ms in range(chunk.shape[0]):
                ms_coords = chunk[i_ms, :, :]
                ms_from_coords = from_coords[i_ms, :, :]

                # Compute distances from all atoms (N) to the closest atom in from_coords (K)
                # ms_coords[:, np.newaxis, :]  -> (N, 1, 3)
                # ms_from_coords[np.newaxis, :, :] -> (1, K, 3)
                # diff will be (N, K, 3) after broadcasting
                diff = ms_coords[:, np.newaxis, :] - ms_from_coords[np.newaxis, :, :]

                # Compute the Euclidean norm along the coordinate axis (axis=2 for (N, K, 3))
                distances = np.linalg.norm(diff, axis=2)  # (N, K)

                min_distances = np.min(distances, axis=1)  # (N,)

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
