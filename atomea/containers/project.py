# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from loguru import logger

from atomea.containers import Container, Energy, Ensemble, Quantum, Time
from atomea.stores import ArrayStore, Store, StoreKind, TableStore


class Project(Container):
    """
    The project is the core interface you have for interacting with data through
    atomea.

    ## Cadence

    All data derived from atomistic calculations are on one of two cadences:
    **ensemble** and **microstate**. Both are derived from the conceptual
    understanding of statistical mechanics.

    > [!NOTE]
    > We do not support ensembles where the number of particles
    > change (i.e., grand canonical).

    Both are defined in the `Cadence` enum in atomea.

    ```python
    from atomea.data import Cadence
    ```

    ### Microstate

    ```python
    Cadence.MICROSTATE
    ```

    A microstate is a single, distinguishable configuration of particles (i.e., atoms)
    where the system's thermodynamic variables are unchanged. In our calculations,
    a microstate could be a:

    - frame of a molecular dynamics simulation trajectory;
    - protein-ligand pose in a docking calculation;
    - transition state of a chemical reaction.

    Data that _could_ change between microstates, such as atomistic coordinates,
    energy, docking score, instantaneous temperature or pressure, dipole moment,
    electronic state, etc., are given a cadence of `MICROSTATE`.

    ### Ensemble

    ```python
    Cadence.ENSEMBLE
    ```

    An ensemble is a collection of microstates where "thermodynamic variables"
    (e.g., Hamiltonian, temperature, number of particles) are constant.
    Changes in any of these variables change the ensemble. This also extends
    to calculation parameters that could&mdash;intentionally or
    not&mdash;change properties of that atomistic system (e.g., force field,
    integration algorithm, docking scoring algorithm, barostat set point, etc.).

    ## Stores

    Dimensionality of data determines how we represent, store, and analyze it.
    We define two working categories: scalars and n-dimensional.

    ### Scalars

    Values that have only one dimension with respect to each microstate are always
    stored in tables using DataFrames in their respective columns.
    This includes energies, thermodynamic variables, calculation parameters, and
    other relevant factors. Data cadence has no influence on storage.

    All data should be stored in a way that assumes multiple ensembles and
    microstates will be present. Each table item must include:

    - `ens_id` (`str`): A unique identification label for an ensemble.
        This can be `"1"`, `"default"`, `"exp3829"`, etc.
    - `run_id` (`str`): An unique, independent run within the same ensemble.
        This often arises when running multiple independent molecular simulation
        trajectories with different random seeds.
    - `micro_id` (`uint`): An index specifying a microstate with some
        relationship to order. This can be a frame in a molecular simulation
        trajectories, docking scores from best to worst, optimization steps, etc.

    ### N-dimensional

    Data with more than one value **must** be stored with arrays with the
    appropriate number of dimensions for multiple values&mdash;even if there is
    only one. Data for all ensemble runs are stored in a single array since
    they are theoretically sampled from the same ensemble.

    Data must also be stored in the same order as the table indices of that
    `Container`. Thus, row indices from tables can be used to slice arrays.
    Note that row indices can change between `Containers` since not all
    data is collected on the same cadence.
    """

    def __init__(
        self,
        array_store: ArrayStore,
        table_store: TableStore,
        label: str = "prj",
    ) -> None:
        """
        Args:
            array_store: Storage backend for all arrays.
            table_store: Storage backend for all tables.
            id: Unique ID for this container.
        """
        assert isinstance(array_store, ArrayStore)
        assert isinstance(table_store, TableStore)
        assert isinstance(label, str)

        self.label = label
        self._stores: dict[StoreKind, Store] = {
            StoreKind.ARRAY: array_store,
            StoreKind.TABLE: table_store,
        }

        self._ensembles: dict[str, Ensemble] = {}

        self.quantum = Quantum(self)

        self.energy = Energy(self)

        self.time = Time(self)

    def __getitem__(self, ens_id: str) -> Ensemble:
        return self.get_ensemble(ens_id)

    def add_ensemble(self, ens_id: str) -> Ensemble:
        """
        Create and register a new Ensemble with given ID, using the project's stores.

        Args:
            ens_id: Unique label for the ensemble.

        Returns:
            The newly created Ensemble.
        """
        assert isinstance(ens_id, str)
        if ens_id in self._ensembles:
            raise ValueError(f"Ensemble '{ens_id}' already exists")

        ens = Ensemble(ens_id=ens_id, parent=self)

        logger.info("Registering ensemble: {}", ens_id)
        self._ensembles[ens_id] = ens
        return ens

    def get_ensemble(self, ens_id: str) -> Ensemble | None:
        """
        Retrieve an existing Ensemble by its ID, or None if not found.
        """
        return self._ensembles.get(ens_id, None)

    def remove_ensemble(self, ens_id: str) -> None:
        """
        Remove an Ensemble from the project by its ID.

        TODO: Need to implement dropping tables and arrays.

        Raises:
            KeyError: if the ensemble does not exist.
        """
        del self._ensembles[ens_id]
        logger.info("Removed ensemble: {}", ens_id)

    def list_ensembles(self) -> list[str]:
        """
        Return all ensemble IDs managed by this project.
        """
        return list(self._ensembles.keys())

    def __repr__(self) -> str:
        return f"<Project ensembles={self.list_ensembles()!r}>"
