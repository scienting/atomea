from typing import Literal

from pydantic import BaseModel, Field

from .amber import AmberSchema


class SimulationSchema(BaseModel):
    compute_platform: Literal["mpi", "cuda"] = Field(default="mpi")
    r"""Options for architecture to run simulations on."""

    amber: AmberSchema = Field(default_factory=AmberSchema)
