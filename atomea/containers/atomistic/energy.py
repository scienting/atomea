import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data, Metadata
from atomea.stores import StoreKind


class Energy(AtomeaContainer):
    """Information that characterizes various energies of a microstates."""

    def __init__(self):
        self._parent = None

    electronic: adt.DataFrame = Data[adt.DataFramePL](
        meta=Metadata(
            uuid="9e4bdf45-0150-4605-9528-e23aed0be9f2",
            cadence=Cadence.MICROSTATE,
            store=StoreKind.TABLE,
            description="Electronic energy",
        ),
        default=None,
    )
    """Electronic energy is the quantum-mechanical expectation value of the electronic Hamiltonian.
    It encompasses all kinetic and potential contributions from the Coulombic interactions of the electrons.
    """

    potential_mm: adt.DataFrame = Data[adt.DataFramePL](
        meta=Metadata(
            uuid="399ff4fb-1b3d-41a4-a87a-8143c1646b28",
            cadence=Cadence.MICROSTATE,
            store=StoreKind.TABLE,
            description="Classical potential energy",
        ),
        default=None,
    )
    """Classical (i.e., molecular mechanics) potential energy including the sum of
    bonded and non-bonded terms.
    """

    kinetic: adt.DataFrame = Data[adt.DataFramePL](
        meta=Metadata(
            uuid="0095592c-587d-4a65-a7f0-d85b588bf2dc",
            cadence=Cadence.MICROSTATE,
            store=StoreKind.TABLE,
            description="Total kinetic energy",
        ),
        default=None,
    )
    """The total kinetic energy of the atomistic system associated with the translational,
    rotational, and vibrational motion of its particles.
    """
