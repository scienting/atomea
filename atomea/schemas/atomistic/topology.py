from typing import Annotated

from pydantic import BaseModel, Field

from ..io import IOBase


class TopologySchema(BaseModel, IOBase):
    """Information that specifies the physical atomistic system."""

    ff_atom_type: Annotated[
        list[str] | None,
        {"cadence": "ensemble", "uuid": "e34c0e1b-0eaa-4679-b060-3fcfe737aa15"},
    ] = Field(default=None)
    """In the context of force fields used in molecular dynamics simulations, an
    "atom type" refers to a specific classification assigned to individual atoms within
    a molecular system based on certain characteristics.
    Atom types play a crucial role in defining the parameters and potential
    energy functions used to calculate forces and motions during a simulation.

    Each atom in a molecular system is assigned a particular atom type, which is
    typically associated with a set of parameters defining its behavior.
    These parameters include values such as atomic mass, partial charges, van der
    Waals radii, and bond, angle, and dihedral force constants.
    The specific values for these parameters are determined based on experimental
    data and quantum mechanical calculations.

    **Cadence:** `ensemble`

    **UUID:** `e34c0e1b-0eaa-4679-b060-3fcfe737aa15`
    """
