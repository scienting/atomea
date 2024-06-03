from typing import Annotated

from pydantic import BaseModel, Field

from ..io import YamlIO


class TopologySchema(BaseModel, YamlIO):
    """Information that specifies the physical atomistic system."""

    ids_entity: Annotated[
        list[int] | None,
        {"cadence": "molecule", "uuid": "b490f2db-548e-4c92-a71a-8222c041ca54"},
    ] = Field(default=None)
    """A uniquely identifying integer specifying what atoms belong to which entities.
    Entities can be a related set of atoms, molecules, or functional group.
    For example, a water and methanol molecule could be `[0, 0, 0, 1, 1, 1, 1, 1, 1]`.

    **Cadence:** `molecule`

    **UUID:** `b490f2db-548e-4c92-a71a-8222c041ca54`
    """

    ids_component: Annotated[
        list[str] | None,
        {"cadence": "molecule", "uuid": "cf39af62-d372-4747-a431-cf2fa0c8e119"},
    ] = Field(default=None)
    """Relates [`ids_entity`][schemas.atomistic.topology.TopologySchema.ids_entity]
    to a fragment label for chemical components or species.
    Labels could be `WAT` or `h2o` for water, `MeOH` for methanol, `bz` for benzene,
    etc. There are no standardized labels for species.

    **Cadence:** `molecule`

    **UUID:** `cf39af62-d372-4747-a431-cf2fa0c8e119`
    """

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
