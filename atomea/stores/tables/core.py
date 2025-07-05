from typing import Any

from abc import ABC, abstractmethod
from pathlib import Path

import polars as pl

from atomea.stores import ArrayDiskFormats, DiskFormat, Store, StoreKind


class TableStore(Store, ABC):
    """
    Abstract interface for tabular data.
    """

    kind = StoreKind.TABLE

    def __init__(
        self,
        path: Path | str,
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
        path: Path | str,
        ensemble_id: str | None = None,
        run_id: int | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
        **kwargs: Any,
    ) -> Any:
        """
        Query a named table using a filter expression.

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

        Returns:
            Filtered DataFrame.
        """
