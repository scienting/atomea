from typing import Any

import sys

import numpy as np
import numpy.typing as npt

try:
    import zarr

    HAS_ZARR = True
except ImportError:
    HAS_ZARR = False
from loguru import logger

from ..manager import StorageManager


class ZarrManager(StorageManager):
    @classmethod
    def get_store(cls, path_file: str) -> str:
        logger.debug(f"Getting zarr store at {path_file}")
        if not HAS_ZARR:
            logger.critical("zarr is not installed!")
            sys.exit(0)
        return path_file

    @classmethod
    def initialize_group(cls, store: Any, key_group: str) -> Any:
        """Check to see if a group at `key_group` exists in `path_file`. If not,
        then it creates a group.

        Args:
            store: A store used directly from the data backend.
            key_group: Key to the group starting from root at `path_file`.
        """
        logger.debug(f"Getting group at {key_group}")
        return zarr.open_group(store=store, path=key_group, mode="a")

    @classmethod
    def append_array(
        cls,
        data: npt.NDArray[Any],
        store: Any,
        key_array: str,
        **kwargs: dict[str, Any],
    ) -> Any:
        logger.debug(f"Appending array at {key_array}")
        parent_group_path = "/".join(key_array.split("/")[:-1])
        parent_group = cls.initialize_group(store, parent_group_path)
        array_name = key_array.split("/")[-1]

        if array_name not in parent_group:
            logger.debug("Did not find the array; creating it now")
            dtype = data.dtype
            if dtype.kind == "U":  # Unicode
                array = parent_group.create_array(
                    name=array_name,
                    shape=data.shape,
                    dtype=data.dtype,
                )
                array[:] = data
            elif dtype.kind == "O":  # Object
                array = parent_group.create_array(
                    name=array_name,
                    shape=data.shape,
                    dtype=np.str_,
                )
                array[:] = data.astype(np.str_)
            else:
                array = parent_group.create_array(
                    name=array_name, shape=data.shape, dtype=data.dtype
                )
                array[:] = data
        else:
            logger.debug("Found the array; appending data")
            parent_group[array_name].append(data)
            array = parent_group[array_name]
        return array
