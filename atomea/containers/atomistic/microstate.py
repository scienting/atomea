import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data, Metadata
from atomea.stores import StoreKind


class Microstates(AtomeaContainer):
    """Information that specifies the physical atomistic microstates."""

    def __init__(self):
        self._parent = None

    atom_z: adt.OptionalUInt8 = Data[adt.UInt8NP](
        meta=Metadata(
            uuid="d051abd9-c815-40b1-ab2d-e7a50a2d3259",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.ARRAY,
            description="Atomic numbers",
        ),
        default=None,
    )
    """The atomic number is a fundamental property of an atom and is denoted by the
    symbol $Z$. It is defined as the number of protons in the nucleus of an atom.
    In a neutral atom, the atomic number also corresponds to the number of electrons
    orbiting the nucleus.
    """

    atom_symbol: adt.OptionalStr = Data[adt.StrNP](
        meta=Metadata(
            uuid="81c21a83-4b72-48c6-a576-4541b468eb90",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.ARRAY,
            description="Elemental symbols",
        ),
        default=None,
    )
    """Elemental symbol based on [`atom_z`]
    [containers.atomistic.microstate.Microstates.atom_z].
    """

    coordinates: adt.OptionalFloat64 = Data[adt.Float64NP](
        meta=Metadata(
            uuid="81c7cec9-beec-4126-b6d8-91bee28951d6",
            cadence=Cadence.MICROSTATE,
            store=StoreKind.ARRAY,
            description="Atomic coordinates",
        ),
        default=None,
    )
    """Coordinates refer to the specific three-dimensional positions of particles
    defined using a set of Cartesian coordinates ($x$, $y$, $z$).
    """

    charge_net: adt.DataFrame = Data[adt.DataFramePL](
        meta=Metadata(
            uuid="6ff82a49-4666-4cbb-978a-409bfa6a511",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Net charge",
        ),
        default=None,
    )
    """The net charge of an atomic system is the overall charge determined by the
    balance between positively charged protons and negatively charged electrons.
    """

    multiplicity: adt.DataFrame = Data[adt.DataFramePL](
        meta=Metadata(
            uuid="8e3eb55a-ed81-46d3-9f34-0ea00fa8c8e4",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Spin multiplicity",
        ),
        default=None,
    )
    """The degeneracy or the number of possible spin states associated with a particular electronic state of a molecule.
    The multiplicity is denoted by the symbol $2S+1$, where $S$ is the total electron spin angular momentum.
    Here, $S$ can take non-negative half-integer values, such as 0, 1/2, 1, 3/2, and so on.
    The multiplicity is always an integer.
    """
