import atomea.typing as adt
from atomea.containers import Container
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Interval(Container):
    """
    A generic interval-based time schema for simulation outputs.
    """

    coord = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="ecd65483-ef07-4893-8b39-5d46118ce97a",
        description="Number of integration steps between writing coordinates",
    )
    """Number of integration steps between writing coordinates.
    """

    energy = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="86ecb388-2760-41af-a989-433402dfcf44",
        description="Number of integration steps between writing energies",
    )
    """Number of integration steps between writing energies.
    """

    velocity = Data[adt.UInt32](
        store_kind=StoreKind.TABLE,
        uuid="2ec3f6c5-01d8-4565-bc50-c91bda41f28c",
        description="Number of integration steps between writing velocities",
    )
    """Number of integration steps between writing velocities.
    """

    def __init__(self, parent: object) -> None:
        self.label = "interval"
        self.cadence = Cadence.RUN
        self._parent = parent

        self.coord.bind_to_container(self)
        self.energy.bind_to_container(self)
        self.velocity.bind_to_container(self)
