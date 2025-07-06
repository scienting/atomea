from typing import TYPE_CHECKING

import numpy as np

from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class AtomTypeIs(SelectionExpression):
    """Selects atoms whose atom type is in the specified list."""

    def __init__(self, atom_types: list[str]):
        self.atom_types = atom_types
        self._cached_all_atom_types: np.ndarray | None = (
            None  # Cache for ensemble-level data
        )

    def evaluate(
        self, ensemble: "Ensemble", run_id: str | None = None
    ) -> np.ndarray:
        # Atom types are ENSEMBLE cadence, so load once and cache
        if self._cached_all_atom_types is None:
            self._cached_all_atom_types = ensemble.topology.atoms.types.values(
                run_id=run_id
            )
            if self._cached_all_atom_types is None:
                coords = ensemble.coordinates.values(
                    run_id=run_id,
                )
                num_atoms = coords.shape[1] if coords is not None else 0
                return np.zeros(num_atoms, dtype=bool)

        return np.isin(self._cached_all_atom_types, self.atom_types)
