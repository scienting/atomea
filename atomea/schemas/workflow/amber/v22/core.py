from pydantic import Field

from atomea.schemas.workflow.amber import AmberSchemaBase
from atomea.schemas.workflow.amber.v22 import (
    Amber22CLI,
    Amber22Forcefield,
    Amber22Inputs,
)


class Amber22Schema(AmberSchemaBase):
    r"""Amber 22 schema for simulation contexts."""

    inputs: Amber22Inputs = Field(default_factory=Amber22Inputs)

    cli: Amber22CLI = Field(default_factory=Amber22CLI)

    ff: Amber22Forcefield = Field(default_factory=Amber22Forcefield)
