from typing import TYPE_CHECKING

import numpy as np

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
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        # Molecule IDs are ENSEMBLE cadence, so load once and cache
        if self._cached_all_mol_ids is None:
            self._cached_all_mol_ids = ensemble.topology.ids.molecules.values(
                ens_id=ensemble.label,
                run_id=run_id,  # Run ID might be relevant for some ensemble data if partitioned
            )
            if self._cached_all_mol_ids is None:
                # Handle case where data might not exist for this ensemble/run
                # Determine number of atoms from coordinates if possible
                coords = ensemble.coordinates.values(
                    run_id=run_id,
                )
                num_atoms = coords.shape[1] if coords is not None else 0
                return np.zeros(num_atoms, dtype=bool)

        return np.isin(self._cached_all_mol_ids, self.mol_ids)
