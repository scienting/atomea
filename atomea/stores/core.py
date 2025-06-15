from typing import Any

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
import polars as pl

from atomea.stores import DiskFormat


class Store(ABC):
    """
    Abstract interface for data stores.
    """

    _store: dict[str, Any]
    disk_format = DiskFormat

    def __init__(self, disk_format: DiskFormat = DiskFormat.NONE) -> None:
        self.disk_format = disk_format

    @abstractmethod
    def write(
        self,
        name: str,
        data: pl.DataFrame | npt.NDArray[np.generic],
    ) -> None:
        """
        Write data.

        Args:
            name: logical name.
            data: Data to write.
        """
        ...

    @abstractmethod
    def append(
        self,
        name: str,
        data: pl.DataFrame | npt.NDArray[np.generic],
    ) -> None:
        """
        Append data.

        Args:
            name: logical name.
            data: Data to append.
        """
        ...

    @abstractmethod
    def read(self, name: str) -> pl.DataFrame | npt.NDArray[np.generic] | None:
        """
        Read the entire data structure with the given name.

        Args:
            name: logical name.

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
