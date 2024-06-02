from typing import Literal

from pydantic import BaseModel, Field

from ...io import YamlIO
from .amber import Amber22Schema


class SimulationSchema(BaseModel, YamlIO):
    compute_platform: Literal["mpi", "cuda"] = Field(default="mpi")
    r"""Options for architecture to run simulations on."""

    amber: Amber22Schema = Field(default_factory=Amber22Schema)
