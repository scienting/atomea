from typing import TYPE_CHECKING, Iterator

import numpy as np
import numpy.typing as npt
from loguru import logger

from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class AtomTypeIs(SelectionExpression):
    """Selects atoms whose atom type is in the specified list."""

    def __init__(self, atom_types: list[str]):
        self.atom_types = atom_types
        self._cached_all_atom_types: npt.NDArray[np.generic] | None = None

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[npt.NDArray[np.bool]]:
        # Atom types are ENSEMBLE cadence, so load once and cache
        if self._cached_all_atom_types is None:
            self._cached_all_atom_types = ensemble.topology.atoms.types.read(
                run_id=run_id
            )

        # If we don't have any atom types, everything should be false.
        if self._cached_all_atom_types is None:
            logger.warning("Did not find any atom types; eval will be false.")
            slice_n_atoms = (0, slice(None), slice(None))
            coords = ensemble.coordinates.read(
                view=slice_n_atoms,
                run_id=run_id,
            )
            num_atoms = coords.shape[1] if coords is not None else 0
            mask = np.zeros(num_atoms, dtype=np.dtype(np.bool))
        else:
            mask = np.isin(self._cached_all_atom_types, self.atom_types)
        yield mask
