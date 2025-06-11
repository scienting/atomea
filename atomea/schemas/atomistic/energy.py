from typing import Annotated

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, field_validator

from atomea.schemas import YamlIO


class EnergySchema(BaseModel, YamlIO):
    """Information that characterizes various energies of a system."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    electronic: Annotated[
        npt.NDArray[np.float64] | None,
        {"cadence": "microstate", "uuid": "9e4bdf45-0150-4605-9528-e23aed0be9f2"},
    ] = Field(default=None)
    """Electronic energy is the quantum-mechanical expectation value of the electronic Hamiltonian.
    It encompasses all kinetic and potential contributions from the Coulombic interactions of the electrons.
    
    **Cadence:** `microstate`

    **UUID:** `9e4bdf45-0150-4605-9528-e23aed0be9f2`
    """

    potential_mm: Annotated[
        npt.NDArray[np.float64] | None,
        {"cadence": "microstate", "uuid": "399ff4fb-1b3d-41a4-a87a-8143c1646b28"},
    ] = Field(default=None)
    """Classical (i.e., molecular mechanics) potential energy including the sum of
    bonded and non-bonded terms.
    
    **Cadence:** `microstate`

    **UUID:** `399ff4fb-1b3d-41a4-a87a-8143c1646b28`
    """

    kinetic: Annotated[
        npt.NDArray[np.float64] | None,
        {"cadence": "microstate", "uuid": "0095592c-587d-4a65-a7f0-d85b588bf2dc"},
    ] = Field(default=None)
    """The total kinetic energy of the atomistic system associated with the translational,
    rotational, and vibrational motion of its particles.
    
    **Cadence:** `microstate`

    **UUID:** `0095592c-587d-4a65-a7f0-d85b588bf2dc`
    """
