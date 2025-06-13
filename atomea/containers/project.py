from loguru import logger

from atomea.containers import AtomeaContainer, Ensemble
from atomea.containers.atomistic import Energy, Quantum, Time
from atomea.stores import ArrayStore, StoreKind, TableStore
from atomea.stores.arrays import NumpyArrayStore
from atomea.stores.tables import PolarsTableStore


class Project(AtomeaContainer):
    """
    Root object managing ensembles and global per-ensemble or per-microstate data.
    """

    def __init__(
        self,
        array_store: ArrayStore | None = None,
        table_store: TableStore | None = None,
    ) -> None:
        """
        Args:
            array_store: Storage backend for all arrays. Defaults to NumpyArrayStore.
            table_store: Storage backend for all tables. Defaults to PolarsTableStore.
        """
        if array_store is None:
            array_store = NumpyArrayStore()
        if table_store is None:
            table_store = PolarsTableStore()
        self._stores: dict[StoreKind, object] = {
            StoreKind.ARRAY: array_store,
            StoreKind.TABLE: table_store,
        }
        self.ensembles: dict[str, Ensemble] = {}

        self.quantum = Quantum()
        self.quantum._parent = self  # type: ignore

        self.energy = Energy()
        self.energy._parent = self  # type: ignore

        self.time = Time()
        self.time._parent = self  # type: ignore

    def add_ensemble(self, ensemble_id: str) -> Ensemble:
        """
        Create and register a new Ensemble with given ID, using the project's stores.

        Args:
            ensemble_id: Unique label for the ensemble.

        Returns:
            The newly created Ensemble.
        """
        assert isinstance(ensemble_id, str)
        if ensemble_id in self.ensembles:
            raise ValueError(f"Ensemble '{ensemble_id}' already exists")

        ens = Ensemble(ensemble_id=ensemble_id, parent=self)  # type: ignore

        logger.info("Registering ensemble: {}", ensemble_id)
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
            logger.info("Removed ensemble: {}", ensemble_id)
        except KeyError:
            raise KeyError(f"Ensemble '{ensemble_id}' not found")

    def list_ensembles(self) -> list[str]:
        """
        Return all ensemble IDs managed by this project.
        """
        return list(self.ensembles.keys())

    def __repr__(self) -> str:
        return f"<Project ensembles={self.list_ensembles()!r}>"
