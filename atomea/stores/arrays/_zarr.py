from typing import Literal

import numpy as np
import numpy.typing as npt
import zarr
from zarr.core.array import Array
from zarr.core.group import Group

from atomea.stores import DiskFormat
from atomea.stores.arrays import ArrayStore


class ZarrArrayStore(ArrayStore):
    """
    Zarr-based array store where each logical path maps to a Zarr array.

    Paths use '/' to denote nested logical structure.
    """

    def __init__(
        self,
        disk_format: DiskFormat,
        store: str | zarr.abc.store.Store,  # type: ignore
        mode: str = "r",
        *args,
        **kwargs,
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
        assert disk_format == DiskFormat.ZARR
        self._store: Group = zarr.open_group(store=store, mode=mode, *args, **kwargs)
        super().__init__(disk_format)

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
        self._store.create_array(
            name=path,
            shape=shape,
            dtype=dtype,
            chunks=chunks,
            overwrite=overwrite * args,
            **kwargs,
        )

    def write(
        self,
        path: str,
        array: npt.NDArray[np.generic],
        slices: tuple[slice, ...] | dict[int, tuple[slice, ...]] | None = None,
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
        z = self._store.get(path=path)
        z.set_basic_selection(slices, array)  # type: ignore

    def append(self, path: str, array: npt.NDArray[np.generic]) -> None:
        """
        Append data along the first axis to an existing Zarr array;
        creates the array with an unlimited first dimension if it does not exist.
        """
        arr_0 = self.read(path)
        arr_0.append(array)  # type: ignore

    def read(
        self,
        path: str,
        slices: tuple[slice, ...] | dict[int, slice | tuple[slice, ...]] | None = None,
    ) -> npt.NDArray[np.generic] | None:
        """
        Read the array from Zarr, optionally returning a subset efficiently.
        Returns None if the path does not exist.
        """
        z = self._store.get(path)
        if slices is None:
            return z[:]  # type: ignore
        if isinstance(slices, dict):
            idx = [slice(None)] * z.ndim  # type: ignore
            for axis, sl in slices.items():
                idx[axis] = sl  # type: ignore
            return z[tuple(idx)]  # type: ignore
        return z[slices]  # type: ignore

    def available(self) -> list[str]:
        """
        List all stored Zarr array paths (joined by '/').
        """
        keys = list(self._store.array_keys())
        return keys
