from atomea.schemas import Ensemble
from atomea.schemas.atomistic import Energy, Quantum, Time
from atomea.schemas.field import StoreKind
from atomea.stores.arrays import ArrayStore
from atomea.stores.tables import TableStore


class Project(Quantum, Energy, Time):
    """
    Root object managing multiple ensembles and global per-ensemble or per-microstate data.

    Inherits schema fields from Quantum, Energy, and Time mixins.
    Ensemble instances also use these stores for per-ensemble and per-microstate data.
    """

    def __init__(
        self,
        arrays: ArrayStore,
        tables: TableStore,
    ) -> None:
        self._stores: dict[StoreKind, object] = {
            StoreKind.ARRAY: arrays,
            StoreKind.TABLE: tables,
        }
        self.ensembles: dict[str, Ensemble] = {}

    def add_ensemble(self, ensemble_id: str) -> Ensemble:
        """
        Create and register a new Ensemble with given ID, using the project's stores.

        Returns:
            The newly created Ensemble.
        """
        if ensemble_id in self.ensembles:
            raise ValueError(f"Ensemble '{ensemble_id}' already exists")
        ens = Ensemble(
            ensemble_id=ensemble_id,
            arrays=self._stores[StoreKind.ARRAY],
            tables=self._stores[StoreKind.TABLE],
        )
        self.ensembles[ensemble_id] = ens
        return ens

    def get_ensemble(self, ensemble_id: str) -> Ensemble | None:
        """
        Retrieve an existing Ensemble by its ID, or None if not found.
        """
        return self.ensembles.get(ensemble_id)

    def remove_ensemble(self, ensemble_id: str) -> None:
        """
        Remove an Ensemble from the project by its ID.

        Raises:
            KeyError: if the ensemble does not exist.
        """
        try:
            del self.ensembles[ensemble_id]
        except KeyError:
            raise KeyError(f"Ensemble '{ensemble_id}' not found")

    def list_ensembles(self) -> list[str]:
        """
        Return all ensemble IDs managed by this project.
        """
        return list(self.ensembles.keys())

    def __repr__(self) -> str:
        return f"<Project ensembles={self.list_ensembles()!r}>"
