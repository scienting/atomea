from enum import Enum


class Cadence(Enum):
    """Frequency at which a field's data changes."""

    ENSEMBLE = 0
    MICROSTATE = 1


class StoreKind(Enum):
    ARRAY = 0
    TABLE = 1


class FieldMeta:
    uuid: str
    cadence: Cadence
    store: StoreKind
    description: str | None

    def __init__(
        self,
        *,
        uuid: str,
        cadence: Cadence,
        store: StoreKind,
        description: str | None = None,
    ) -> None:
        self.uuid = uuid
        self.cadence = cadence
        self.store = store
        self.description = description
