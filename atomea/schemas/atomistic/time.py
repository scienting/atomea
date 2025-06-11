from typing import Annotated

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, field_validator

from atomea.schemas import YamlIO


class TimeSchema(BaseModel, YamlIO):
    """
    A generic interval-based time schema for simulation outputs.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    time_step: Annotated[
        float | None,
        {"cadence": "ensemble", "uuid": "ec042cd8-c4de-4655-b663-cb96493b2ded"},
    ] = Field(default=None)
    """Integration time step in femtoseconds (fs).
    
    **Cadence:** `ensemble`

    **UUID:** `ec042cd8-c4de-4655-b663-cb96493b2ded`
    """

    interval_coord: Annotated[
        int | None,
        {"cadence": "ensemble", "uuid": "ecd65483-ef07-4893-8b39-5d46118ce97a"},
    ] = Field(default=None)
    """Number of integration steps between writing coordinates.
    
    **Cadence:** `ensemble`

    **UUID:** `ecd65483-ef07-4893-8b39-5d46118ce97a`
    """

    interval_energy: Annotated[
        int | None,
        {"cadence": "ensemble", "uuid": "86ecb388-2760-41af-a989-433402dfcf44"},
    ] = Field(default=None)
    """Number of integration steps between writing energies.
    
    **Cadence:** `ensemble`

    **UUID:** `86ecb388-2760-41af-a989-433402dfcf44`
    """

    interval_velocity: Annotated[
        int | None,
        {"cadence": "ensemble", "uuid": "2ec3f6c5-01d8-4565-bc50-c91bda41f28c"},
    ] = Field(default=None)
    """Number of integration steps between writing velocities.
    
    **Cadence:** `ensemble`

    **UUID:** `2ec3f6c5-01d8-4565-bc50-c91bda41f28c`
    """
