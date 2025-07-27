# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

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
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        self._store: dict[str, pl.DataFrame] = {}
        assert disk_format not in ArrayDiskFormats or disk_format == DiskFormat.NONE
        super().__init__(path, mode=mode, disk_format=disk_format, **kwargs)

    @abstractmethod
    def query(
        self,
        path: Path | str,
        ens_id: str | None = None,
        run_id: str | None = None,
        micro_id: int | None = None,
        filter_expr: str | None = None,
        **kwargs: Any,
    ) -> Any:
        """
        Query a named table using a filter expression.

        Args:
            path: Container/table name.
            ens_id: A unique identification label for an ensemble.
                This can be `"1"`, `"default"`, `"exp3829"`, etc.
            run_id: An unique, independent run within the same ensemble.
                This often arises when running multiple independent molecular
                simulation trajectories with different random seeds.
            micro_id: An index specifying a microstate with some relationship to
                order. This can be a frame in a molecular simulation trajectories,
                docking scores from best to worst, optimization steps, etc.
            filter_expr: string expression to filter rows.

        Returns:
            Filtered DataFrame.
        """
