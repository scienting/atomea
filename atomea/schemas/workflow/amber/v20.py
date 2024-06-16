from typing import Literal

from pydantic import Field

from ..ff import ForcefieldSchemaBase
from .cli import AmberCLIBase
from .inputs import AmberInputsBase
from .schema import AmberSchemaBase


class Amber20Inputs(AmberInputsBase):
    pass


class Amber20CLI(AmberCLIBase):
    pass


class Amber20Forcefield(ForcefieldSchemaBase):
    protein: (
        Literal[
            "ff19SB",
            "ff14SB",
            "ff99SB",
            "ff15ipq",
            "fb15",
            "ff03ua",
        ]
        | None
    ) = Field(default=None)
    r"""Options for protein force fields.

    -   [ff19SB](https://md.crumblearn.org/mm/examples/protein/sb/19/)
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


class Amber20Schema(AmberSchemaBase):
    r"""Amber 20 schema for simulation contexts."""

    inputs: Amber20Inputs = Field(default_factory=Amber20Inputs)

    cli: Amber20CLI = Field(default_factory=Amber20CLI)

    ff: Amber20Forcefield = Field(default_factory=Amber20Forcefield)
