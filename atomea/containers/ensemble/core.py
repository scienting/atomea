# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import TYPE_CHECKING

import atomea.typing as adt
from atomea.containers import Container
from atomea.containers.ensemble import Topology
from atomea.data import Cadence, Data, OptionalSliceSpec
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

    def __init__(self, ens_id: str, parent: "Project") -> None:
        self.label = ens_id
        self.cadence = Cadence.MICROSTATE
        self._parent = parent

        self.topology = Topology(self)

        self.coordinates.bind_to_container(self)

        self._n_micro: int | None = None

    def __repr__(self) -> str:
        return f"<Ensemble id={self.label!r}>"

    def n_micro(self, view: OptionalSliceSpec = None, run_id: str | None = None):
        """Total number of microstates."""
        if self._n_micro:
            return self._n_micro

        to_check = (self.coordinates,)
        for data_check in to_check:
            if data_check.store_kind != StoreKind.ARRAY:
                continue
            data = data_check.read(view=view, run_id=run_id)
            if data is not None:
                self._n_micro = int(data.shape[0])
                return self._n_micro

        return None
