from typing import Any

from abc import ABC
from pathlib import Path

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
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        assert disk_format in ArrayDiskFormats or disk_format == DiskFormat.NONE
        super().__init__(path, disk_format=disk_format, **kwargs)
