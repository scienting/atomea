from typing import Any

import os
from pathlib import Path

import polars as pl
from loguru import logger

from atomea.data import OptionalSliceSpec
from atomea.stores import DiskFormat
from atomea.stores.tables import TableStore


class PolarsTableStore(TableStore):
    """
    Each table is a Polars DataFrame.
    """

    def __init__(
        self,
        path: Path | str,
        disk_format: DiskFormat = DiskFormat.PARQUET,
        mode: str = "r",
        **kwargs: Any,
    ) -> None:
        super().__init__(path, disk_format=disk_format)

    @classmethod
    def check_columns(cls, columns):
        if (
            "ensemble_id" not in columns
            or "run_id" not in columns
            or "microstate_id" not in columns
        ):
            raise ValueError(
                "Table must include 'ensemble_id', 'run_id', 'microstate_id' columns"
            )

    def write(
        self,
        path: Path | str,
        data: pl.DataFrame,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        """
        Write a Polars DataFrame to the named table.

        The DataFrame must contain `ensemble_id`, `run_id`, and `microstate_id` columns.
        """
        self.check_columns(data.columns)
        self._store[str(path)] = data

    def append(self, path: Path | str, data: pl.DataFrame, **kwargs: Any) -> None:
        """
        Append a Polars DataFrame to the named table.

        The DataFrame must contain `ensemble_id` and `microstate_id` columns.
        """
        path = str(path)
        self.check_columns(data.columns)
        self._store[path] = pl.concat([self._store[path], data], how="vertical")

    def get(
        self,
        path: Path | str,
        **kwargs: Any,
    ) -> pl.DataFrame | None:
        """Get the store-specific object that represents the data stored here.

        For example, if the data is stored on disk using some type of memory map, this
        would return the memory map object, not the in-memory data. If you want to
        guarantee the data is loaded into memory, use `read`.
        """
        try:
            return self._store[str(path)]
        except KeyError:
            logger.warning(f"Table '{path}' not found! Returning empty DataFrame")
            return pl.DataFrame()

    def read(
        self, path: Path | str, view: OptionalSliceSpec = None, **kwargs: Any
    ) -> pl.Series:
        """
        Read the entire table with the given name.

        TODO: Consider returning only the series of the exact column instead of the
        whole dataframe.
        """
        df = self.get(path)
        return df

    def query(
        self,
        path: Path | str,
        ensemble_id: str | None = None,
        run_id: int | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
        **kwargs: Any,
    ) -> pl.DataFrame:
        """
        Query a named table by keys and/or additional expression.

        Args:
            path: Container/table name.
            ensemble_id: A unique identification label for an ensemble.
                This can be `"1"`, `"default"`, `"exp3829"`, etc.
            run_id: An unique, independent run within the same ensemble.
                This often arises when running multiple independent molecular
                simulation trajectories with different random seeds.
            microstate_id: An index specifying a microstate with some relationship to
                order. This can be a frame in a molecular simulation trajectories,
                docking scores from best to worst, optimization steps, etc.
            filter_expr: string expression to filter rows.
        """
        df = self.get(path)

        # Early return if df is empty
        if df.shape == (0, 0):
            return df

        if ensemble_id is not None:
            df = df.filter(pl.col("ensemble_id") == ensemble_id)
        if run_id is not None:
            df = df.filter(pl.col("run_id") == run_id)
        if microstate_id is not None:
            df = df.filter(pl.col("microstate_id") == microstate_id)
        if filter_expr:
            df = df.filter(pl.parse_expr(filter_expr))
        return df

    def available(self) -> list[str]:
        """List all table names."""
        return list(self._store.keys())

    def dump(self, **kwargs: Any) -> None:
        """
        Dump all stored tables to files in the specified directory/prefix.
        """
        for name, df in self._store.items():
            path = os.path.join(self.path, f"{name}")
            if self.disk_format == DiskFormat.CSV:
                path += ".csv"
                df.write_csv(path, **kwargs)
            elif self.disk_format == DiskFormat.PARQUET:
                path += ".parquet"
                df.write_parquet(path, **kwargs)
            elif self.disk_format == DiskFormat.IPC:
                path += ".ipc"
                df.write_ipc(path, **kwargs)
            elif self.disk_format == DiskFormat.EXCEL:
                path += ".xlsx"
                df.write_excel(path, **kwargs)
            else:
                raise ValueError(f"Unsupported file_type: {self.disk_format}")
