from typing import TYPE_CHECKING

import atomea.typing as adt
from atomea.containers import Container
from atomea.containers.ensemble import Topology
from atomea.data import Cadence, Data
from atomea.stores import StoreKind

if TYPE_CHECKING:
    from atomea.containers import Project


class Ensemble(Container):
    """
    The `Ensemble` class represents a collection of molecular structures,
    each referred to as a microstate. This class is used to
    manage and validate an ensemble of molecular data, facilitating the handling of
    multiple molecular configurations, such as those produced during atomistic
    calculations.

    Only data that could reasonably change shape or dimensions between ensembles
    (due to different numbering or ordering of atoms) should be stored here. All other
    data should be stored in a [`Project`][schemas.Project].
    """

    coordinates = Data[adt.Float64](
        store_kind=StoreKind.ARRAY,
        uuid="81c7cec9-beec-4126-b6d8-91bee28951d6",
        description="Atomic coordinates",
    )
    """Coordinates refer to the specific three-dimensional positions of particles
    defined using a set of Cartesian coordinates ($x$, $y$, $z$).
    """

    def __init__(self, ensemble_id: str, parent: "Project") -> None:
        self.label = ensemble_id
        self.cadence = Cadence.MICROSTATE
        self._parent = parent

        self.topology = Topology(self)

        self.coordinates.bind_to_container(self)

    def __repr__(self) -> str:
        return f"<Ensemble id={self.label!r}>"
