from pydantic import BaseModel, Field

from .id import IdentificationSchema
from .qc import QCSchema
from .system import SystemSchema


class MoleculeSchema(BaseModel):
    identification: IdentificationSchema = Field(default_factory=IdentificationSchema)

    qc: QCSchema = Field(default_factory=QCSchema)

    system: SystemSchema = Field(default_factory=SystemSchema)


class EnsembleSchema(BaseModel):
    frames: list[MoleculeSchema] = Field(default=[])
