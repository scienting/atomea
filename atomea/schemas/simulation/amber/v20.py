"""Simulation contexts for Amber"""
from pydantic import Field

from .inputs import AmberInputsBase
from .schema import AmberSchemaBase


class Amber20Inputs(AmberInputsBase):
    pass


class Amber20Schema(AmberSchemaBase):
    r"""Amber 20 schema for simulation contexts."""

    inputs: Amber20Inputs = Field(default_factory=Amber20Inputs)
