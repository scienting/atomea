from pydantic import BaseModel, Field

from ..io import YamlIO
from .molecule import MoleculeSchema


class EnsembleSchema(BaseModel, YamlIO):
    """The EnsembleSchema class is a Pydantic model designed to represent a
    collection of molecular structures, referred to as frames. This class is used to
    manage and validate an ensemble of molecular data, facilitating the handling of
    multiple molecular configurations, such as those produced during molecular dynamics
    simulations or geometry optimizations.
    """

    frames: list[MoleculeSchema] = Field(default=[])
    """The frames attribute is a list that stores instances of MoleculeSchema. Each
    instance represents a single molecular structure or configuration within the
    ensemble. This attribute allows the EnsembleSchema to manage multiple molecular
    structures collectively, making it easier to handle data from simulations that
    produce multiple frames, such as molecular dynamics trajectories or conformational
    scans."""
