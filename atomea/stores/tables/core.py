from typing import Any

import os
from abc import ABC, abstractmethod

import polars as pl

from atomea.stores import ArrayDiskFormats, DiskFormat, Store


class TableStore(Store, ABC):
    """
    Abstract interface for tabular data mappable to ensembles and microstates,
    e.g., energy, time stamps, indices, and other per-microstate properties.
    """

    def __init__(
        self,
        path: str,
        *args: Any,
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            store: Directory to store all tables. We recommend "<project>.tables".
            disk_format: Format to store data on disk.
        """
        self._store: dict[str, pl.DataFrame] = {}
        assert disk_format not in ArrayDiskFormats
        super().__init__(path, *args, disk_format=disk_format, **kwargs)

    @abstractmethod
    def query(
        self,
        name: str,
        ensemble_id: str | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
        *kwargs: Any,
    ) -> pl.DataFrame:
        """
        Query a named table using a filter expression.

        Args:
            name: logical table name.
            filter_expr: string expression to filter rows, e.g., "ensemble_id=='e1'".

        Returns:
            Filtered DataFrame.
        """
        ...
