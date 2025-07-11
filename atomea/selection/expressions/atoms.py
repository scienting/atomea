from typing import TYPE_CHECKING, Iterator

import numpy as np
import numpy.typing as npt

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class SelectByAtomType(SelectionExpression):
    """Selects atoms whose atom type is in the specified list."""

    def __init__(self, atom_types: list[str]):
        self.atom_types = atom_types
        self._cached_all_atom_types: npt.NDArray[np.generic] | None = None

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
        if self._cached_all_atom_types is None:
            self._cached_all_atom_types = ensemble.topology.atoms.types.read(
                run_id=run_id
            )

        if self._cached_all_atom_types is None:
            mask: adt.Bool = self.init_mask()
        else:
            mask: adt.Bool = np.isin(self._cached_all_atom_types, self.atom_types)
        yield mask
