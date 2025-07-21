# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import TYPE_CHECKING

from collections.abc import Iterator

import numpy as np

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.selection.expressions import SelectionExpression

if TYPE_CHECKING:
    from atomea.containers import Ensemble


class SelectByMoleculeID(SelectionExpression):
    """Selects atoms whose molecule ID is in the specified list."""

    def __init__(self, mol_ids: list[int]):
        self.mol_ids = mol_ids
        self._cached_all_mol_ids: adt.UInt32 | None = None

    def evaluate(
        self,
        ensemble: "Ensemble",
        run_id: str | None = None,
        micro_id: OptionalSliceSpec = None,
    ) -> Iterator[adt.Bool]:
        if self._cached_all_mol_ids is None:
            self._cached_all_mol_ids = ensemble.topology.ids.molecules.read(
                ens_id=ensemble.label,
                run_id=run_id,
            )

        if self._cached_all_mol_ids is None:
            mask: adt.Bool = self.init_mask()
        else:
            mask: adt.Bool = np.isin(self._cached_all_mol_ids, self.mol_ids)
        yield mask
