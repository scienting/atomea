from typing import Literal

from pydantic import Field

from .cli import AmberCLIBase
from .inputs import AmberInputsBase
from .schema import AmberSchemaBase


class Amber22Inputs(AmberInputsBase):
    pass


class Amber22CLI(AmberCLIBase):
    pass


class Amber22Schema(AmberSchemaBase):
    r"""Amber 22 schema for simulation contexts."""

    inputs: Amber22Inputs = Field(default_factory=Amber22Inputs)

    ff_protein: Literal[
        "ff19SB",
        "ff14SB",
        "ff99SB",
        "ff15ipq",
        "fb15",
        "ff03ua",
    ] = Field(default="ff19SB")
    r"""Options for protein force fields.

    ### ff19SB

    DOI:  [10.1021/acs.jctc.9b00591](https://doi.org/10.1021/acs.jctc.9b00591)

    ### ff14SB

    DOI:  [10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)

    ### ff99SB

    DOI:  [10.1002/prot.21123](https://doi.org/10.1002/prot.21123)

    ### ff15ipq

    DOI:  [10.1021/acs.jctc.6b00567](https://doi.org/10.1021/acs.jctc.6b00567)

    ### fb15

    DOI:  [10.1021/acs.jpcb.7b02320](https://doi.org/10.1021/acs.jpcb.7b02320)

    ### ff03ua

    DOI:  [10.1021/jp060163v](https://doi.org/10.1021/jp060163v)
    """

    ff_water: Literal[
        "tip4p",
        "tip4pew",
        "tip5p",
        "spce",
        "spceb",
        "opc",
        "opc3",
        "opc3pol",
        "opc3pol",
        "pol3",
        "tip3pfb",
        "tip4pfb",
    ] = Field(default="opc3")
    r"""Options for water force fields."""
