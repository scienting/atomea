import atomea.typing as adt
from atomea.containers import Container
from atomea.data import Cadence, Data, OptionalSliceSpec
from atomea.stores import StoreKind


class Bonds(Container):
    """Information that specifies (constant) binding."""

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
        self.label = "bonds"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.covalent.bind_to_container(self)

        self._n_covalent_bonds = None

    def n_covalent_bonds(
        self, view: OptionalSliceSpec = None, run_id: str | None = None
    ):
        """Total number of covalent bonds."""
        if self._n_covalent_bonds:
            return self._n_covalent_bonds

        to_check = (self.covalent,)
        for data_check in to_check:
            if data_check.store_kind != StoreKind.ARRAY:
                continue
            data = data_check.read(view=view, run_id=run_id)
            if data is not None:
                self._n_covalent_bonds = int(data.shape[0])
                return self._n_covalent_bonds

        return None
