from typing import Literal

from pydantic import Field

from ..ff import ForcefieldSchemaBase
from .cli import AmberCLIBase
from .inputs import AmberInputsBase
from .schema import AmberSchemaBase


class Amber18Inputs(AmberInputsBase):
    pass


class Amber18CLI(AmberCLIBase):
    pass


class Amber18Forcefield(ForcefieldSchemaBase):
    protein: (
        Literal[
            "ff14SB",
            "ff99SB",
            "ff15ipq",
            "fb15",
            "ff03ua",
        ]
        | None
    ) = Field(default=None)
    r"""Options for protein force fields.

    -   [ff14SB](https://md.crumblearn.org/mm/examples/protein/sb/14/)
    -   [ff99SB](https://md.crumblearn.org/mm/examples/protein/sb/99/)
    -   [ff15ipq](https://md.crumblearn.org/mm/examples/protein/ipq/15/)
    -   [fb15](https://md.crumblearn.org/mm/examples/protein/fb/15/)
    -   ff03ua
    """

    water: (
        Literal[
            "tip4p",
            "tip4pew",
            "tip5p",
            "spce",
            "spceb",
            "opc",
            "opc3",
            "pol3",
            "tip3pfb",
            "tip4pfb",
        ]
        | None
    ) = Field(default=None)
    r"""Options for water force fields."""


class Amber18Schema(AmberSchemaBase):
    r"""Amber 18 schema for simulation contexts."""

    inputs: Amber18Inputs = Field(default_factory=Amber18Inputs)

    cli: Amber18CLI = Field(default_factory=Amber18CLI)

    ff: Amber18Forcefield = Field(default_factory=Amber18Forcefield)
