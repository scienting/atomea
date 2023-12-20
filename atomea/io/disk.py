from typing import Any

import os
from collections import defaultdict
from collections.abc import Iterable

import numpy as np

from ..schema import Atomea
from .backends._parquet import tabular_parquet
from .backends._zarr import array_zarr


class DiskData:
    """Manage data with a single backend."""

    def __init__(self, backend_array=array_zarr, backend_tabular=tabular_parquet):
        self.backend_array = backend_array
        self.backend_tabular = backend_tabular

    def store(
        self, dest: str, data: dict[str, Any] | Iterable[dict[str, Any]], atomea: Atomea
    ) -> None:
        """Store data into array and tabular formats into a single location. Only
        includes keys that are shared between all `data`"""
        if isinstance(data, dict):
            data = [data]
        common_keys = set(data[0].keys())
        for d in data:
            common_keys = common_keys.intersection(d.keys())

        merged_dict = defaultdict(list)
        for d in data:
            for key in common_keys:
                merged_dict[key].append(d[key])
        merged_dict = dict(merged_dict)

        with atomea as schema:
            data_tabular = {}
            for k, v in merged_dict.items():
                if isinstance(v[0], (np.ndarray,)):
                    v = np.array(v)
                    data_path = os.path.join(dest, k)

                    self.backend_array(data_path, v, dtype=schema[k]["dtype"])
                else:
                    data_tabular[k] = v
            data_path = os.path.join(dest, "tabular")
            self.backend_tabular(data_path, data_tabular)
