import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Atoms(AtomeaContainer):
    """Information that specifies atoms."""

    atomic_numbers = Data[adt.UInt8](
        store_kind=StoreKind.ARRAY,
        uuid="d051abd9-c815-40b1-ab2d-e7a50a2d3259",
        description="Atomic numbers",
    )
    """The atomic number is a fundamental property of an atom and is denoted by the
    symbol $Z$. It is defined as the number of protons in the nucleus of an atom.
    In a neutral atom, the atomic number also corresponds to the number of electrons
    orbiting the nucleus.
    """

    symbols = Data[adt.Str](
        store_kind=StoreKind.ARRAY,
        uuid="81c21a83-4b72-48c6-a576-4541b468eb90",
        description="Elemental symbols",
    )
    """Elemental symbol based on [`atom_z`]
    [containers.atomistic.microstate.Microstates.atom_z].
    """

    types = Data[adt.Str](
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
        self.label = "atoms"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent

        self.atomic_numbers.bind_to_container(self)
        self.symbols.bind_to_container(self)
        self.types.bind_to_container(self)
