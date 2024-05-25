import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, field_validator


class SystemSchema(BaseModel):
    """Information that specifies the physical atomistic system."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    atom_z: npt.NDArray[np.uint8] | None = Field(default=None)
    """The atomic number is a fundamental property of an atom and is denoted by the
    symbol $Z$. It is defined as the number of protons in the nucleus of an atom.
    In a neutral atom, the atomic number also corresponds to the number of electrons
    orbiting the nucleus."""

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

    atom_symbol: list[str] | None = Field(default=None)
    """The atomic number is a fundamental property of an atom and is denoted by the
    symbol $Z$. It is defined as the number of protons in the nucleus of an atom.
    In a neutral atom, the atomic number also corresponds to the number of electrons
    orbiting the nucleus."""

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

    coordinates: npt.NDArray[np.float64 | np.float32] | None = Field(default=None)
    """Coordinates refer to the specific three-dimensional positions of particles
    defined using a set of Cartesian coordinates ($x$, $y$, $z$).
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

    charge_net: int | None = Field(default=None)
    """The net charge of an atomic system is the overall charge determined by the
    balance between positively charged protons and negatively charged electrons."""

    multiplicity: int | None = Field(default=None)
    """The degeneracy or the number of possible spin states associated with a particular electronic state of a molecule.
    The multiplicity is denoted by the symbol $2S+1$, where $S$ is the total electron spin angular momentum.
    Here, $S$ can take non-negative half-integer values, such as 0, 1/2, 1, 3/2, and so on.
    The multiplicity is always an integer"""

    ids_entity: list[int] | None = Field(default=None)
    """A uniquely identifying integer specifying what atoms belong to which entities.
    Entities can be a related set of atoms, molecules, or functional group.
    For example, a water and methanol molecule could be `[0, 0, 0, 1, 1, 1, 1, 1, 1]`.
    """

    ids_component: list[str] | None = Field(default=None)
    """Relates ``entity_id`` to a fragment label for chemical components or species.
    Labels could be WAT or h2o for water, MeOH for methanol, bz for benzene, etc.
    There are no standardized labels for species."""
