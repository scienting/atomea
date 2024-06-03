from typing import Literal

from pydantic import BaseModel, Field

from ...io import YamlIO


class SimulationSchema(BaseModel, YamlIO):
    compute_platform: Literal["mpi", "cuda"] = Field(default="mpi")
    r"""Options for architecture to run simulations on."""
