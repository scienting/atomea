"""Simulation contexts for Amber"""
from pydantic import Field

from .inputs import AmberInputsBase
from .schema import AmberSchemaBase


class Amber18Inputs(AmberInputsBase):
    pass


class Amber18Schema(AmberSchemaBase):
    r"""Amber 18 schema for simulation contexts."""

    inputs: Amber18Inputs = Field(default_factory=Amber18Inputs)
