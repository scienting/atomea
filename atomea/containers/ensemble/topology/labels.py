# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import atomea.typing as adt
from atomea.containers import Container
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Labels(Container):
    """Information that specifies the physical atomistic ensemble.

    Note that this topology is assumed constant; a reactive topology needs
    to be implemented.
    """

    components = Data[adt.Str](
        store_kind=StoreKind.ARRAY,
        uuid="3466ffde-1ac0-4e07-ad5f-832420c3943f",
        description="Maps a Component ID to a human-readable label.",
    )
    """An array of human-readable component labels on a per-atom basis."""

    def __init__(self, parent: object) -> None:
        self.label = "labels"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.components.bind_to_container(self)
