# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import atomea.typing as adt
from atomea.containers import Container
from atomea.containers.ensemble.topology.connectivity import Bonds
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Connectivity(Container):
    """Information that specifies (constant) binding."""

    angles = Data[adt.UInt64](
        store_kind=StoreKind.ARRAY,
        uuid="876fff34-27c8-4caf-8fe9-94bb1c149ead",
        description="Atom indices of angles with respect to covalent bonds",
    )
    """
    Specifies atom indices that form angles from covalent bonds.
    """

    dihedrals = Data[adt.UInt64](
        store_kind=StoreKind.ARRAY,
        uuid="3bbd4368-869e-463e-86b0-3d9833f89255",
        description="Atom indices of dihedrals with respect to covalent bonds",
    )
    """
    Specifies atom indices that form dihedrals from covalent bonds.
    """

    def __init__(self, parent: object) -> None:
        self.label = "bonds"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.bonds = Bonds(self)

        self.angles.bind_to_container(self)
        self.dihedrals.bind_to_container(self)
