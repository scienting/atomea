import numpy as np
import numpy.typing as npt

from atomea.data import Cadence, Data, Metadata
from atomea.stores import StoreKind


class Microstates:
    """Information that specifies the physical atomistic microstates."""

    def __init__(self):
        self._parent = None

    atom_z: npt.NDArray[np.uint8] | None = Data[npt.NDArray[np.uint8]](
        dtype=np.uint8,
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

    **Cadence:** `ENSEMBLE`

    **UUID:** `d051abd9-c815-40b1-ab2d-e7a50a2d3259`
    """

    atom_symbol: npt.NDArray[np.str_] | None = Data[npt.NDArray[np.str_]](
        dtype=str,
        meta=Metadata(
            uuid="81c21a83-4b72-48c6-a576-4541b468eb90",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.ARRAY,
            description="Elemental symbols",
        ),
        default=None,
    )
    """Elemental symbol based on [`atom_z`]
    [schemas.atomistic.microstates.Microstates.atom_z].

    **Cadence:** `ENSEMBLE`

    **UUID:** `81c21a83-4b72-48c6-a576-4541b468eb90`
    """

    coordinates: npt.NDArray[np.float64] | None = Data[npt.NDArray[np.float64]](
        dtype=float,
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

    **Cadence:** `MICROSTATE`

    **UUID:** `81c7cec9-beec-4126-b6d8-91bee28951d6`
    """

    charge_net: npt.NDArray[np.int8] | None = Data[npt.NDArray[np.int8]](
        dtype=np.int8,
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

    **Cadence:** `ENSEMBLE`

    **UUID:** `6ff82a49-4666-4cbb-978a-409bfa6a511`
    """

    multiplicity: npt.NDArray[np.uint8] | None = Data[npt.NDArray[np.uint8]](
        dtype=np.uint8,
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

    **Cadence:** `ENSEMBLE`

    **UUID:** `8e3eb55a-ed81-46d3-9f34-0ea00fa8c8e4`
    """
