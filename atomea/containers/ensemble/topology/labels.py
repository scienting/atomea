import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Labels(AtomeaContainer):
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
