from typing import Any, Literal

import os

import polars as pl

from atomea.stores.tables import TableStore


class PolarsTableStore(TableStore):
    """
    In-memory store for narrow tables of microstate data keyed by separate columns:
    `ensemble_id` and `microstate_id`.

    Each table is a Polars DataFrame with at least these two key columns.

    TODO: Need to handle reading tables from disk; don't need to keep everything in memory.
    """

    def write(self, name: str, table: pl.DataFrame, append: bool = False) -> None:
        """
        Write or append a Polars DataFrame to the named table.

        The DataFrame must contain `ensemble_id` and `microstate_id` columns.
        """
        if "ensemble_id" not in table.columns or "microstate_id" not in table.columns:
            raise ValueError(
                "Table must include 'ensemble_id' and 'microstate_id' columns"
            )
        if append and name in self._store:
            self._store[name] = pl.concat([self._store[name], table], how="vertical")
        else:
            self._store[name] = table

    def read(self, name: str) -> pl.DataFrame:
        """
        Read the entire table with the given name.

        Raises KeyError if the table does not exist.
        """
        try:
            return self._store[name]
        except KeyError:
            raise KeyError(f"Table '{name}' not found")

    def query(
        self,
        name: str,
        ensemble_id: str | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
    ) -> pl.DataFrame:
        """
        Query a named table by keys and/or additional expression.

        - If `ensemble_id` or `microstate_id` is provided, filters on those.
        - `filter_expr` is a Polars expression string for further filtering.
        """
        df = self.read(name)
        if ensemble_id is not None:
            df = df.filter(pl.col("ensemble_id") == ensemble_id)
        if microstate_id is not None:
            df = df.filter(pl.col("microstate_id") == microstate_id)
        if filter_expr:
            df = df.filter(pl.parse_expr(filter_expr))
        return df

    def available(self) -> list[str]:
        """
        List all table names.
        """
        return list(self._store.keys())

    def dump(
        self,
        prefix: str,
        file_type: Literal["csv", "parquet", "ipc", "excel"],
        **options: Any,
    ) -> None:
        """
        Dump all stored tables to files in the specified directory/prefix.

        Args:
            prefix: directory or file prefix where files will be written.
            file_type: one of 'csv', 'parquet', or 'ipc' (Arrow IPC).
            options: backend-specific write options (e.g., compression='snappy').
        """
        os.makedirs(prefix, exist_ok=True)
        for name, df in self._store.items():
            path = os.path.join(prefix, f"{name}.{file_type}")
            if file_type == "csv":
                df.write_csv(path, **options)
            elif file_type == "parquet":
                df.write_parquet(path, **options)
            elif file_type == "ipc":
                df.write_ipc(path, **options)
            elif file_type == "excel":
                df.write_excel(path, **options)
            else:
                raise ValueError(f"Unsupported file_type: {file_type}")
