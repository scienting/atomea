from typing import TYPE_CHECKING, Iterator

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class MolIdIs(SelectionExpression):
    """Selects atoms whose molecule ID is in the specified list."""

    def __init__(self, mol_ids: list[int]):
        self.mol_ids = mol_ids
        self._cached_all_mol_ids: np.ndarray | None = (
            None  # Cache for ensemble-level data
        )

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        # Molecule IDs are ENSEMBLE cadence, so load once and cache
        if self._cached_all_mol_ids is None:
            self._cached_all_mol_ids = ensemble.topology.ids.molecules.read(
                ens_id=ensemble.label,
                run_id=run_id,  # Run ID might be relevant for some ensemble data if partitioned
            )
            if self._cached_all_mol_ids is None:
                # Handle case where data might not exist for this ensemble/run
                # Determine number of atoms from coordinates if possible
                coords = ensemble.coordinates.read(
                    run_id=run_id,
                )
                num_atoms = coords.shape[1] if coords is not None else 0
                yield np.zeros(num_atoms, dtype=np.dtype(np.bool))

        yield np.isin(self._cached_all_mol_ids, self.mol_ids)
