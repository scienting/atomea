from typing import Any

import os
from abc import ABC, abstractmethod

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.stores import DiskFormat


class Store(ABC):
    """
    Abstract interface for data stores.
    """

    _store: Any

    def __init__(
        self,
        path: str,
        *args: Any,
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.disk_format = disk_format

    @abstractmethod
    def write(
        self,
        path: str,
        data: adt.PassableData,
        *args: Any,
        slice: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        """
        Write data.

        Args:
            path: logical name.
            data: Data to write.
        """
        ...

    @abstractmethod
    def append(
        self,
        path: str,
        data: adt.PassableData,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Append data.

        Args:
            path: logical name.
            data: Data to append.
        """
        ...

    @abstractmethod
    def read(
        self, path: str, *args, slice: OptionalSliceSpec = None, **kwargs: Any
    ) -> adt.OptionalPassableData:
        """
        Read the entire data structure with the given name.

        Args:
            path: logical name.

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

    def dump(self, *args, prefix: str = "", **kwargs: Any) -> None:
        """Save all data to disk. Overwrite this method to implement when needed."""
        pass
