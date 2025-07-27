# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import atomea.typing as adt
from atomea.containers import Container
from atomea.containers.time import Interval
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Time(Container):
    """
    A generic interval-based time schema for simulation outputs.
    """

    time_step = Data[adt.UInt8](
        store_kind=StoreKind.TABLE,
        uuid="ec042cd8-c4de-4655-b663-cb96493b2ded",
        description="Integration time step in femtoseconds (fs)",
    )
    """Integration time step in femtoseconds (fs).
    """

    def __init__(self, parent: object) -> None:
        self.label = "time"
        self.cadence = Cadence.RUN
        self._parent = parent

        self.interval = Interval(self)

        self.time_step.bind_to_container(self)
