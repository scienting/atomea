from abc import ABC, abstractmethod

import polars as pl

from atomea.stores import Store, DiskFormat, ArrayDiskFormats


class TableStore(Store, ABC):
    """
    Abstract interface for tabular data mappable to ensembles and microstates,
    e.g., energy, time stamps, indices, and other per-microstate properties.
    """

    def __init__(self, disk_format: DiskFormat = DiskFormat.NONE) -> None:
        self._store: dict[str, pl.DataFrame] = {}
        assert disk_format not in ArrayDiskFormats
        super().__init__(disk_format)

    @abstractmethod
    def write(self, name: str, table: pl.DataFrame, append: bool = False) -> None:
        """
        Write a pandas DataFrame to a named table. Optionally append.

        Args:
            name: logical table name (e.g., "microstate_index").
            table: DataFrame containing columns [ensemble_id, microstate_id, ...].
            append: if True, append rows to an existing table; otherwise overwrite.
        """
        ...

    @abstractmethod
    def read(self, name: str) -> pl.DataFrame | None:
        """
        Read the entire table with the given name.

        Args:
            name: logical table name.

        Returns:
            DataFrame previously stored.
        """
        ...

    @abstractmethod
    def query(
        self,
        name: str,
        ensemble_id: str | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
    ) -> pl.DataFrame:
        """
        Query a named table using a filter expression (e.g., pandas query syntax).

        Args:
            name: logical table name.
            filter_expr: string expression to filter rows, e.g., "ensemble_id=='e1'".

        Returns:
            Filtered DataFrame.
        """
        ...

    @abstractmethod
    def available(self) -> list[str]:
        """
        List all table names currently stored.

        Returns:
            A list of table names.
        """
        ...
