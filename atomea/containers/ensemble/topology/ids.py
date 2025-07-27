# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import atomea.typing as adt
from atomea.containers import Container
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class IDs(Container):
    """Information that specifies the physical atomistic ensemble.

    Note that this topology is assumed constant; a reactive topology needs
    to be implemented.
    """

    molecules = Data[adt.UInt32](
        store_kind=StoreKind.ARRAY,
        uuid="b490f2db-548e-4c92-a71a-8222c041ca54",
        description="Uniquely identifying integer mapping atoms to a molecule.",
    )
    """A uniquely identifying integer specifying which atoms belong to a single,
    physical, covalently bonded molecule. A molecule can be an organic compound,
    ion, protein, single-stranded DNA, etc.
    For example, a water and methanol molecule could be `[0, 0, 0, 1, 1, 1, 1, 1, 1]`.
    """

    components = Data[adt.UInt32](
        store_kind=StoreKind.ARRAY,
        uuid="cf39af62-d372-4747-a431-cf2fa0c8e119",
        description="Unique integer that maps atoms to substructures within a molecule.",
    )
    """Assigns atoms from a single `id_molecule` into various substructures.
    Substructures and represent functional groups in organic compounds,
    amino acids in proteins, nucleotides in DNA, etc.

    The number of unique components is constant between microstates.
    Atoms can change components; for example, a proton transfer would result in
    one hydrogen changing from one component to another (e.g., `32` -> `59`).
    """

    def __init__(self, parent: object) -> None:
        self.label = "ids"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.molecules.bind_to_container(self)
        self.components.bind_to_container(self)
