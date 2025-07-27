# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

from abc import ABC, abstractmethod
from collections.abc import Iterator
from pathlib import Path

from atomea.data import OptionalSliceSpec
from atomea.stores import ArrayDiskFormats, DiskFormat, Store, StoreKind


class ArrayStore(Store, ABC):
    """
    Abstract interface for storing and retrieving arrays,
    e.g., coordinates, velocities, forces.
    """

    kind = StoreKind.ARRAY

    def __init__(
        self,
        path: Path | str,
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        assert disk_format in ArrayDiskFormats or disk_format == DiskFormat.NONE
        super().__init__(path, mode=mode, disk_format=disk_format, **kwargs)

    def shape(self, path: Path | str) -> tuple[int, ...]:
        """Get the shape of the underlying array."""
        z = self.get(path)
        if z is None:
            return None
        return z.shape

    @abstractmethod
    def iter(
        self,
        path: Path | str,
        elements: OptionalSliceSpec = None,
        view: OptionalSliceSpec = None,
        chunk_size: int = 1,
        **kwargs: Any,
    ) -> Iterator[Any]:
        """Yield chunks of data instead of reading all into memory.

        Args:
            path:
            view: View for all but the first dimension.
            chunk: Number of data points of the first axis to yield.
        """
