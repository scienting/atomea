from typing import Any

import os
from abc import ABC, abstractmethod
from pathlib import Path

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.stores import DiskFormat, StoreKind


class Store(ABC):
    """
    An abstract interface for handling data stored in memory or on disk.
    """

    kind: StoreKind
    _store: Any

    def __init__(
        self,
        path: Path | str,
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            path: Path where data will be stored.
            disk_format: File format when writing data to disk.
        """
        os.makedirs(path, exist_ok=True)
        self.path = path
        self.disk_format = disk_format

    @abstractmethod
    def write(
        self,
        path: Path | str,
        data: Any,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        """
        Write data.

        Args:
            path: Path where data will be stored.
            data: Data to write.
        """
        ...

    @abstractmethod
    def append(
        self,
        path: Path | str,
        data: Any,
        **kwargs: Any,
    ) -> None:
        """
        Append data.

        Args:
            path: Path where data will be stored.
            data: Data to append.
        """
        ...

    @abstractmethod
    def read(
        self, path: Path | str, view: OptionalSliceSpec = None, **kwargs: Any
    ) -> adt.OptionalPassableData:
        """
        Read the entire data structure with the given name.

        Args:
            path: Path where data will be stored.

        Returns:
            The requested data.
        """
        ...

    @abstractmethod
    def available(self) -> list[str]:
        """
        List all data names currently stored.

        Returns:
            A list of data names.
        """
        ...

    def dump(self, **kwargs: Any) -> None:
        """Save all data to disk. Overwrite this method to implement when needed."""
        pass
