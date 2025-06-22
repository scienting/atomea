import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Time(AtomeaContainer):
    """
    A generic interval-based time schema for simulation outputs.
    """

    time_step = Data[adt.UInt8](
        store_kind=StoreKind.TABLE,
        uuid="ec042cd8-c4de-4655-b663-cb96493b2ded",
        description="Integration time step in femtoseconds (fs)",
    )
    """Integration time step in femtoseconds (fs).
    """

    interval_coord = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="ecd65483-ef07-4893-8b39-5d46118ce97a",
        description="Number of integration steps between writing coordinates",
    )
    """Number of integration steps between writing coordinates.
    """

    interval_energy = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="86ecb388-2760-41af-a989-433402dfcf44",
        description="Number of integration steps between writing energies",
    )
    """Number of integration steps between writing energies.
    """

    interval_velocity = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="2ec3f6c5-01d8-4565-bc50-c91bda41f28c",
        description="Number of integration steps between writing velocities",
    )
    """Number of integration steps between writing velocities.
    """

    def __init__(self, parent: object) -> None:
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent
        self.time_step.bind_to_container(self)
        self.interval_coord.bind_to_container(self)
        self.interval_energy.bind_to_container(self)
        self.interval_velocity.bind_to_container(self)
