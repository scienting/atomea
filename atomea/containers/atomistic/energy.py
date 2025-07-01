import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Energy(AtomeaContainer):
    """Information that characterizes various energies of a microstates."""

    electronic = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="9e4bdf45-0150-4605-9528-e23aed0be9f2",
        description="Electronic energy",
    )
    """Electronic energy is the quantum-mechanical expectation value of the electronic Hamiltonian.
    It encompasses all kinetic and potential contributions from the Coulombic interactions of the electrons.
    """

    potential_mm = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="399ff4fb-1b3d-41a4-a87a-8143c1646b28",
        description="Classical potential energy",
    )
    """Classical (i.e., molecular mechanics) potential energy including the sum of
    bonded and non-bonded terms.
    """

    kinetic = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="0095592c-587d-4a65-a7f0-d85b588bf2dc",
        description="Total kinetic energy",
    )
    """The total kinetic energy of the atomistic system associated with the translational,
    rotational, and vibrational motion of its particles.
    """

    def __init__(self, parent: object) -> None:
        self.id = "energy"
        self.cadence = Cadence.MICROSTATE
        self._parent = parent
        self.electronic.bind_to_container(self)
        self.potential_mm.bind_to_container(self)
        self.kinetic.bind_to_container(self)
