from typing import Any

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
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            path: Directory to store all tables. We recommend "<project>.tables".
            disk_format: Format to store data on disk.
        """
        self._store: dict[str, pl.DataFrame] = {}
        assert disk_format not in ArrayDiskFormats
        super().__init__(path, disk_format=disk_format, **kwargs)

    @abstractmethod
    def query(
        self,
        path: str,
        ensemble_id: str | None = None,
        run_id: str | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
        **kwargs: Any,
    ) -> Any:
        """
        Query a named table using a filter expression.

        Args:
            path: logical table name.
            filter_expr: string expression to filter rows, e.g., "ensemble_id=='e1'".

        Returns:
            Filtered DataFrame.
        """
        ...
