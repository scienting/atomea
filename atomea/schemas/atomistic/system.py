from typing import Annotated

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, field_validator

from ..io import YamlIO


class SystemSchema(BaseModel, YamlIO):
    """Information that specifies the physical atomistic system."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    atom_z: Annotated[
        npt.NDArray[np.uint8] | None,
        {"cadence": "ensemble", "uuid": "d051abd9-c815-40b1-ab2d-e7a50a2d3259"},
    ] = Field(default=None)
    """The atomic number is a fundamental property of an atom and is denoted by the
    symbol $Z$. It is defined as the number of protons in the nucleus of an atom.
    In a neutral atom, the atomic number also corresponds to the number of electrons
    orbiting the nucleus.

    **Cadence:** `ensemble`

    **UUID:** `d051abd9-c815-40b1-ab2d-e7a50a2d3259`
    """

    @field_validator("atom_z")
    def validate_atom_z(cls, v):
        if v is None:
            return v

        # Convert to numpy array
        if not isinstance(v, np.ndarray):
            v = np.array(v)

        # Validate shape
        if v.ndim != 1:
            raise ValueError("atom_z must have one dimension")

        return v

    atom_symbol: Annotated[
        list[str] | None,
        {"cadence": "ensemble", "uuid": "81c21a83-4b72-48c6-a576-4541b468eb90"},
    ] = Field(default=None)
    """Elemental symbol based on [`atom_z`]
    [schemas.atomistic.system.SystemSchema.atom_z].

    **Cadence:** `ensemble`

    **UUID:** `81c21a83-4b72-48c6-a576-4541b468eb90`
    """

    @field_validator("atom_symbol")
    def validate_atom_symbol(cls, v):
        if v is None:
            return v

        # Convert to numpy array
        if not isinstance(v, np.ndarray):
            v = np.array(v)

        # Validate shape
        if v.ndim != 1:
            raise ValueError("atom_z must have one dimension")

        return v

    coordinates: Annotated[
        npt.NDArray[np.float64 | np.float32] | None,
        {"cadence": "molecule", "uuid": "81c7cec9-beec-4126-b6d8-91bee28951d6"},
    ] = Field(default=None)
    """Coordinates refer to the specific three-dimensional positions of particles
    defined using a set of Cartesian coordinates ($x$, $y$, $z$).

    **Cadence:** `molecule`

    **UUID:** `81c7cec9-beec-4126-b6d8-91bee28951d6`
    """

    @field_validator("coordinates")
    def validate_coordinates(cls, v):
        if v is None:
            return v

        # Convert to numpy array
        if not isinstance(v, np.ndarray):
            v = np.array(v)

        # Validate shape
        if v.ndim != 2 or v.shape[-1] != 3:
            raise ValueError("Coordinates must have shape (n_atoms, 3)")

        return v

    charge_net: Annotated[
        int | None,
        {"cadence": "molecule", "uuid": "6ff82a49-4666-4cbb-978a-409bfa6a511"},
    ] = Field(default=None)
    """The net charge of an atomic system is the overall charge determined by the
    balance between positively charged protons and negatively charged electrons.

    **Cadence:** `molecule`

    **UUID:** `6ff82a49-4666-4cbb-978a-409bfa6a511`
    """

    multiplicity: Annotated[
        int | None,
        {"cadence": "molecule", "uuid": "8e3eb55a-ed81-46d3-9f34-0ea00fa8c8e4"},
    ] = Field(default=None)
    """The degeneracy or the number of possible spin states associated with a particular electronic state of a molecule.
    The multiplicity is denoted by the symbol $2S+1$, where $S$ is the total electron spin angular momentum.
    Here, $S$ can take non-negative half-integer values, such as 0, 1/2, 1, 3/2, and so on.
    The multiplicity is always an integer.

    **Cadence:** `molecule`

    **UUID:** `8e3eb55a-ed81-46d3-9f34-0ea00fa8c8e4`
    """
