import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Topology(AtomeaContainer):
    """Information that specifies the physical atomistic microstates."""

    ids_molecule = Data[adt.UInt32](
        store_kind=StoreKind.ARRAY,
        uuid="b490f2db-548e-4c92-a71a-8222c041ca54",
        description="Uniquely identifying integer mapping atoms to a molecule.",
    )
    """A uniquely identifying integer specifying which atoms belong to a single,
    physical, covalently bonded molecule. A molecule can be an organic compound,
    ion, protein, single-stranded DNA, etc.
    For example, a water and methanol molecule could be `[0, 0, 0, 1, 1, 1, 1, 1, 1]`.
    """

    ids_component = Data[adt.UInt32](
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

    labels_component = Data[adt.Str](
        store_kind=StoreKind.ARRAY,
        uuid="3466ffde-1ac0-4e07-ad5f-832420c3943f",
        description="Maps a Component ID to a human-readable label.",
    )
    """An array of human-readable component labels on a per-atom basis."""

    ff_atom_type = Data[adt.Str](
        store_kind=StoreKind.ARRAY,
        uuid="e34c0e1b-0eaa-4679-b060-3fcfe737aa15",
        description="Classical force field atom type",
    )
    """In the context of force fields used in molecular dynamics simulations, an
    "atom type" refers to a specific classification assigned to individual atoms within
    a molecular system based on certain characteristics.
    Atom types play a crucial role in defining the parameters and potential
    energy functions used to calculate forces and motions during a simulation.

    Each atom in a molecular system is assigned a particular atom type, which is
    typically associated with a set of parameters defining its behavior.
    These parameters include values such as atomic mass, partial charges, van der
    Waals radii, and bond, angle, and dihedral force constants.
    The specific values for these parameters are determined based on experimental
    data and quantum mechanical calculations.
    """

    def __init__(self, parent: object) -> None:
        self.id = "topology"
        self.cadence = Cadence.MICROSTATE
        self._parent = parent
        self.ids_molecule.bind_to_container(self)
        self.ids_component.bind_to_container(self)
        self.labels_component.bind_to_container(self)
        self.ff_atom_type.bind_to_container(self)
