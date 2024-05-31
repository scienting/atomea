from typing import Any

import numpy.typing as npt
import zarr
from loguru import logger
from numcodecs import MsgPack

from ..manager import StorageManager


class ZarrManager(StorageManager):
    @classmethod
    def get_store(cls, path_file: str) -> zarr.DirectoryStore:
        logger.debug(f"Getting zarr store at {path_file}")
        return zarr.DirectoryStore(path_file)

    @classmethod
    def initialize_group(cls, store: Any, key_group: str) -> Any:
        """Check to see if a group at `key_group` exists in `path_file`. If not,
        then it creates a group.

        Args:
            store: A store used directly from the data backend.
            key_group: Key to the group starting from root at `path_file`.
        """
        logger.debug(f"Getting group at {key_group}")
        return zarr.group(store=store).require_group(key_group)

    @classmethod
    def append_array(
        cls, data: npt.NDArray[Any], storage: zarr.DirectoryStore, key_array: str
    ) -> zarr.Array:
        logger.debug(f"Appending array at {key_array}")
        parent_group_path = "/".join(key_array.split("/")[:-1])
        parent_group = cls.initialize_group(storage, parent_group_path)
        array_name = key_array.split("/")[-1]

        if array_name not in parent_group:
            logger.debug("Did not find the array; creating it now")
            dtype = data.dtype
            if dtype.kind in {"U", "O"}:  # Unicode or Object dtype
                array = parent_group.create_dataset(
                    name=array_name,
                    data=data,
                    shape=data.shape,
                    dtype=data.dtype,
                    object_codec=MsgPack(),
                )
            else:
                array = parent_group.create_dataset(
                    name=array_name, data=data, shape=data.shape, dtype=data.dtype
                )
        else:
            logger.debug("Found the array; appending data")
            parent_group[array_name].append(data)
            array = parent_group[array_name]
        return array
