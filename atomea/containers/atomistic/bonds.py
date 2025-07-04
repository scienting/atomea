import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Bonding(AtomeaContainer):
    """Information that specifies (constant) binding.
    """

    covalent = Data[adt.UInt64](
        store_kind=StoreKind.ARRAY,
        uuid="f1e38cb3-2461-4aa6-bc69-a7555074aeb4",
        description="Atom indices of covalent bonds",
    )
    """
    Specifies atom indices that are covalently bonded. Bonds are specified in a two-dimensional
    array where columns represent the atom index involved in the covalent bond.
    """

    def __init__(self, parent: object) -> None:
        self.id = "bonding"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.covalent.bind_to_container(self)
