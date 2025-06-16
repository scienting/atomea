from enum import Enum


class Cadence(Enum):
    """Frequency at which a field's data changes."""

    ENSEMBLE = 0
    MICROSTATE = 1
