from typing import Literal

import numpy as np
import numpy.typing as npt
import zarr
from zarr.core.array import Array
from zarr.core.group import Group

from atomea.stores.arrays import ArrayStore


class ZarrArrayStore(ArrayStore):
    """
    Zarr-based array store where each logical path maps to a Zarr array.

    Paths use '/' to denote nested logical structure.
    """

    def __init__(
        self, store: str | zarr.abc.store.Store, mode: str = "r", *args, **kwargs
    ) -> None:
        """
        Open a Zarr store using
        [`zarr.open_group`](https://zarr.readthedocs.io/en/stable/api/zarr/index.html#zarr.open_group).
        Please refer to it's
        [documentation](https://zarr.readthedocs.io/en/stable/api/zarr/index.html#zarr.open_group)
        for information on available `args` and `kwargs`.

        Args:
            store: Store or path to directory in file system or name of zip file.
                Strings are interpreted as paths on the local file system and used as
                the root argument to `zarr.storage.LocalStore`.
                Dictionaries are used as the store_dict argument in
                `zarr.storage.MemoryStore`.
                By default (store=None) a new zarr.storage.MemoryStore is created.
            mode: Persistence modes:

                - `r` means read only (must exist);
                - `r+` means read/write (must exist);
                - `a` means read/write (create if doesn't exist);
                - `w` means create (overwrite if exists);
                - `w-` means create (fail if exists).
        """
        self._root: Group = zarr.open_group(store=store, mode=mode, *args, **kwargs)

    def create(
        self,
        path: str,
        shape: tuple[int, ...],
        dtype: npt.DTypeLike | None = None,
        chunks: tuple[int, ...] | Literal["auto"] = "auto",
        overwrite: bool = False,
        *args,
        **kwargs,
    ) -> None:
        """
        Pre-allocate an array with the given shape and dtype.

        Args:
            path: hierarchical key, e.g. 'coords'.
            shape: full shape of the array.
            dtype: numpy-compatible dtype.
            chunks: chunk shape or 'auto'.
            overwrite: if True, delete existing before create; otherwise error if exists.
        """
        group_path, name = path.rsplit("/", 1) if "/" in path else ("", path)
        group = self._root.require_group(group_path)
        if name in group:
            if overwrite:
                del group[name]
            else:
                raise ValueError(f"Array '{path}' already exists")
        group.create_array(
            store=self._root,
            name=name,
            shape=shape,
            dtype=dtype,
            chunks=chunks,
            *args,
            **kwargs,
        )

    @staticmethod
    def _write_full(name, group, array, overwrite):
        if name in group:
            if overwrite:
                del group[name]
            else:
                # assume shape matches, else user should call create
                ds: Array = group[name]
                ds[:] = array
                return
        # create new array with same shape
        ds = group.create_array(
            name=name,
            shape=array.shape,
            dtype=array.dtype,
            chunks="auto",
            overwrite=True,
        )
        ds[:] = array

    @staticmethod
    def _write_slice(name, group, array, slices, overwrite):
        ds: Array = group[name]
        if isinstance(slices, dict):
            idx = [slice(None)] * ds.ndim
            for axis, sl in slices.items():
                idx[axis] = sl
            ds[tuple(idx)] = array
        else:
            ds[slices] = array

    def write(
        self,
        path: str,
        array: npt.NDArray[np.generic],
        slices: tuple[slice, ...] | dict[int, tuple[slice, ...]] | None = None,
        overwrite: bool = False,
    ) -> None:
        """
        Write data to an array, whole or sliced.

        Args:
            path: hierarchical key of the array.
            array: data to write.
            slices: if None, overwrite entire array (requires create or existing),
                    else write into the specified slice region.
            overwrite: if True and slices is None, replaces existing array definition.
        """
        group_path, name = path.rsplit("/", 1) if "/" in path else ("", path)
        group = self._root.require_group(group_path)
        if slices is None:
            self._write_full(name, group, array, overwrite)
        else:
            self._write_slice(name, group, array, slices, overwrite)

    def append(self, path: str, array: npt.NDArray[np.generic]) -> None:
        """
        Append data along the first axis to an existing Zarr array;
        creates the array with an unlimited first dimension if it does not exist.
        """
        group_path, name = path.rsplit("/", 1)
        group = self._root.require_group(group_path)
        arr = np.asarray(array)
        if name in group:
            ds: Array = group[name]  # existing Array
            ds.append(arr)
        else:
            ds = group.create_array(
                name=name,
                shape=(0, *arr.shape[1:]),
                dtype=arr.dtype,
                chunk_shape=None,
                overwrite=True,
            )
            ds.append(arr)

    def read(
        self,
        path: str,
        slices: tuple[slice, ...] | dict[int, slice | tuple[slice, ...]] | None = None,
    ) -> npt.NDArray[np.generic] | None:
        """
        Read the array from Zarr, optionally returning a subset efficiently.
        Returns None if the path does not exist.
        """
        group_path, name = path.rsplit("/", 1)
        try:
            group = self._root[group_path]
        except KeyError:
            return None
        if name not in group:
            return None
        ds: Array = group[name]
        if slices is None:
            return ds[:]
        if isinstance(slices, dict):
            idx = [slice(None)] * ds.ndim
            for axis, sl in slices.items():
                idx[axis] = sl
            return ds[tuple(idx)]
        return ds[slices]

    def available(self) -> list[str]:
        """
        List all stored Zarr array paths (joined by '/').
        """
        paths: list[str] = []

        def _recurse(group: Group, prefix: str) -> None:
            for key in group.array_keys():
                subpath = f"{prefix}/{key}" if prefix else key
                paths.append(subpath)
            for sub in group.group_keys():
                subgrp = group[sub]
                subpath = f"{prefix}/{sub}" if prefix else sub
                _recurse(subgrp, subpath)

        _recurse(self._root, "")
        return paths
