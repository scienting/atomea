# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

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
