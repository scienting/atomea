from loguru import logger

from atomea.containers import AtomeaContainer, Ensemble
from atomea.containers.atomistic import Energy, Quantum, Time
from atomea.stores import ArrayStore, Store, StoreKind, TableStore


class Project(AtomeaContainer):
    """
    Root object managing ensembles and global per-ensemble or per-microstate data.
    """

    def __init__(
        self,
        array_store: ArrayStore,
        table_store: TableStore,
    ) -> None:
        """
        Args:
            array_store: Storage backend for all arrays.
            table_store: Storage backend for all tables.
        """
        assert isinstance(array_store, ArrayStore)
        assert isinstance(table_store, TableStore)

        self._stores: dict[StoreKind, Store] = {
            StoreKind.ARRAY: array_store,
            StoreKind.TABLE: table_store,
        }
        """
        We store array and table stores as a dictionary mainly to use `_stores`
        as an attribute only existing in a Project class.
        """
        self.ensembles: dict[str, Ensemble] = {}

        self.quantum = Quantum(self)

        self.energy = Energy(self)

        self.time = Time(self)

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

        ens = Ensemble(ensemble_id=ensemble_id, parent=self)

        logger.info("Registering ensemble: {}", ensemble_id)
        self.ensembles[ensemble_id] = ens
        return ens

    def get_ensemble(self, ensemble_id: str) -> Ensemble | None:
        """
        Retrieve an existing Ensemble by its ID, or None if not found.
        """
        return self.ensembles.get(ensemble_id, None)

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
