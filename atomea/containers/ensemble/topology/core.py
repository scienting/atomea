from atomea.containers import Container
from atomea.containers.ensemble.topology import Atoms, Connectivity, IDs, Labels
from atomea.data import Cadence


class Topology(Container):
    """Information that specifies the physical atomistic ensemble.

    Note that this topology is assumed constant; a reactive topology needs
    to be implemented.
    """

    def __init__(self, parent: object) -> None:
        self.label = "topology"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.atoms = Atoms(self)
        self.connectivity = Connectivity(self)
        self.ids = IDs(self)
        self.labels = Labels(self)
