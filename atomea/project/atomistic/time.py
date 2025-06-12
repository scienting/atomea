import numpy as np
import numpy.typing as npt

from atomea.data import Cadence, Data, Metadata
from atomea.stores import StoreKind


class Time:
    """
    A generic interval-based time schema for simulation outputs.
    """

    def __init__(self):
        self._parent = None

    time_step: Data[npt.NDArray[np.uint8]] = Data[npt.NDArray[np.uint8]](
        dtype=np.uint8,
        meta=Metadata(
            uuid="ec042cd8-c4de-4655-b663-cb96493b2ded",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Integration time step in femtoseconds (fs)",
        ),
        default=None,
    )
    """Integration time step in femtoseconds (fs).

    **Cadence:** `ENSEMBLE`

    **UUID:** `ec042cd8-c4de-4655-b663-cb96493b2ded`
    """

    interval_coord: npt.NDArray[np.uint32] | None = Data[npt.NDArray[np.uint32]](
        dtype=np.uint32,
        meta=Metadata(
            uuid="ecd65483-ef07-4893-8b39-5d46118ce97a",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Number of integration steps between writing coordinates",
        ),
        default=None,
    )
    """Number of integration steps between writing coordinates.

    **Cadence:** `ENSEMBLE`

    **UUID:** `ecd65483-ef07-4893-8b39-5d46118ce97a`
    """

    interval_energy: npt.NDArray[np.uint32] | None = Data[npt.NDArray[np.uint32]](
        dtype=np.uint32,
        meta=Metadata(
            uuid="86ecb388-2760-41af-a989-433402dfcf44",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Number of integration steps between writing energies",
        ),
        default=None,
    )
    """Number of integration steps between writing energies.

    **Cadence:** `ENSEMBLE`

    **UUID:** `86ecb388-2760-41af-a989-433402dfcf44`
    """

    interval_velocity: npt.NDArray[np.uint32] | None = Data[npt.NDArray[np.uint32]](
        dtype=np.uint32,
        meta=Metadata(
            uuid="2ec3f6c5-01d8-4565-bc50-c91bda41f28c",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Number of integration steps between writing velocities",
        ),
        default=None,
    )
    """Number of integration steps between writing velocities.

    **Cadence:** `ENSEMBLE`

    **UUID:** `2ec3f6c5-01d8-4565-bc50-c91bda41f28c`
    """
