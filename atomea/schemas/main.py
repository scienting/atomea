from pydantic import BaseModel, Field

from .atomistic import EnsembleSchema
from .id import IdentificationSchema
from .slurm import SlurmSchema


class AtomeaSchema(BaseModel):
    """Base class for all Atomea schemas."""

    slurm: SlurmSchema = Field(default_factory=SlurmSchema)
    identification: IdentificationSchema = Field(default_factory=IdentificationSchema)
    ensemble: EnsembleSchema = Field(default_factory=EnsembleSchema)
