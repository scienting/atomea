from atomea.schemas.atomistic import Microstates, Topology
from atomea.schemas.field import StoreKind
from atomea.stores.arrays import ArrayStore
from atomea.stores.tables import TableStore


class Ensemble(Microstates, Topology):
    """
    The `Ensemble` class represents a collection of molecular structures,
    each referred to as a microstate. This class is used to
    manage and validate an ensemble of molecular data, facilitating the handling of
    multiple molecular configurations, such as those produced during atomistic
    calculations.

    Only data that could reasonably change shape or dimensions between ensembles
    (due to different numbering or ordering of atoms) should be stored here. All other
    data should be stored in a [`Project`][schemas.Project].
    """

    def __init__(
        self, ensemble_id: str, arrays: ArrayStore, tables: TableStore
    ) -> None:
        self.id: str = ensemble_id

        self._stores: dict[StoreKind, object] = {
            StoreKind.ARRAY: arrays,
            StoreKind.TABLE: tables,
        }

    def __repr__(self) -> str:
        return f"<Ensemble id={self.id!r}>"
