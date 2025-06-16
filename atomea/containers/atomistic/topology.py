import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Topology(AtomeaContainer):
    """Information that specifies the physical atomistic microstates."""

    ids_entity = Data[adt.UInt32](
        cadence=Cadence.MICROSTATE,
        store_kind=StoreKind.ARRAY,
        uuid="b490f2db-548e-4c92-a71a-8222c041ca54",
        description="Uniquely identifying integer mapping atoms to chemical entities",
    )
    """A uniquely identifying integer specifying what atoms belong to which entities.
    Entities can be a related set of atoms, molecules, or functional group.
    For example, a water and methanol molecule could be `[0, 0, 0, 1, 1, 1, 1, 1, 1]`.
    """

    ids_component = Data[adt.Str](
        cadence=Cadence.ENSEMBLE,
        store_kind=StoreKind.ARRAY,
        uuid="cf39af62-d372-4747-a431-cf2fa0c8e119",
        description="Fragments label for chemical entities",
    )
    """Relates [`ids_entity`][schemas.atomistic.topology.Topology.ids_entity]
    to a fragment label for chemical components or species.
    Labels could be `WAT` or `h2o` for water, `MeOH` for methanol, `bz` for benzene,
    etc. There are no standardized labels for species.
    """

    ff_atom_type = Data[adt.Str](
        cadence=Cadence.ENSEMBLE,
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
        self._parent = parent
        self.ids_entity.bind_to_container(self)
        self.ids_component.bind_to_container(self)
        self.ff_atom_type.bind_to_container(self)
