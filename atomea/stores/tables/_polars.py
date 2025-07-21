# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

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
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.PARQUET,
        **kwargs: Any,
    ) -> None:
        super().__init__(path, mode=mode, disk_format=disk_format)

    @classmethod
    def check_columns(cls, columns):
        if (
            "ens_id" not in columns
            or "run_id" not in columns
            or "micro_id" not in columns
        ):
            raise ValueError(
                "Table must include 'ens_id', 'run_id', 'micro_id' columns"
            )

    def write(
        self,
        path: Path | str,
        data: pl.DataFrame,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        if self.mode == "r":
            raise ValueError("Cannot create when in 'r' mode")
        self.check_columns(data.columns)
        self._store[str(path)] = data

    def append(self, path: Path | str, data: pl.DataFrame, **kwargs: Any) -> None:
        if self.mode in ("r", "w-"):
            raise ValueError(f"Cannot append when in '{self.mode}' mode")
        path = str(path)
        self.check_columns(data.columns)
        self._store[path] = pl.concat([self._store[path], data], how="vertical")

    def get(
        self,
        path: Path | str,
        **kwargs: Any,
    ) -> pl.DataFrame:
        try:
            return self._store[str(path)]
        except KeyError:
            logger.warning(f"Table '{path}' not found! Returning empty DataFrame")
            return pl.DataFrame()

    def read(
        self, path: Path | str, view: OptionalSliceSpec = None, **kwargs: Any
    ) -> pl.Series:
        df = self.get(path)
        return df

    def query(
        self,
        path: Path | str,
        ens_id: str | None = None,
        run_id: int | None = None,
        micro_id: int | None = None,
        filter_expr: str | None = None,
        **kwargs: Any,
    ) -> pl.DataFrame:
        df = self.get(path)

        # Early return if df is empty
        if df.shape == (0, 0):
            return df

        if ens_id is not None:
            df = df.filter(pl.col("ens_id") == ens_id)
        if run_id is not None:
            df = df.filter(pl.col("run_id") == run_id)
        if micro_id is not None:
            df = df.filter(pl.col("micro_id") == micro_id)
        if filter_expr:
            df = df.filter(pl.parse_expr(filter_expr))
        return df

    def available(self) -> list[str]:
        return list(self._store.keys())

    def flush(self, **kwargs: Any) -> None:
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
