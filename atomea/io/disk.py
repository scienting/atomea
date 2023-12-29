from typing import Any

import os
from collections import defaultdict
from collections.abc import Iterable

import numpy as np

from ..schema import Atomea
from ..utils import get_obj_from_string


class DiskData:
    """Manage data with a single backend."""

    def __init__(self, backend_array="zarr", backend_tabular="parquet"):
        self.backend_array = get_obj_from_string("atomea.io.backends._" + backend_array)
        self.backend_tabular = get_obj_from_string(
            "atomea.io.backends._" + backend_tabular
        )

    def store(
        self, dest: str, data: dict[str, Any] | Iterable[dict[str, Any]], atomea: Atomea
    ) -> None:
        """Store data into array and tabular formats into a single location. Only
        includes keys that are shared between all `data`"""
        if isinstance(data, dict):
            data = [data]
        common_keys = set(data[0].keys())  # type: ignore
        for d in data:
            common_keys = common_keys.intersection(d.keys())

        merged_dict = defaultdict(list)
        for d in data:
            for key in common_keys:
                merged_dict[key].extend(d[key])
        merged_dict = dict(merged_dict)  # type: ignore

        with atomea as schema:
            data_tabular: dict[str, dict[str, Any]] = {}
            for k, v in merged_dict.items():
                if isinstance(v[0], np.ndarray) and not schema[k]["tabular"]:
                    v = np.array(v)  # type: ignore
                    data_path = os.path.join(dest, k)

                    self.backend_array.array_write(
                        data_path, v, dtype=schema[k]["dtype"]
                    )
                else:
                    length: str = schema[k]["length"]
                    if length not in data_tabular.keys():
                        data_tabular[length] = {}
                    data_tabular[length][k] = v

            if len(data_tabular) > 0:
                for length_key, tab_data in data_tabular.items():
                    data_path = os.path.join(dest, length_key)
                    self.backend_tabular.tabular_write(
                        data_path, tab_data, atomea.fields(tab_data.keys())
                    )
