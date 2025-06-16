from typing import Any, Literal

import numpy as np
import numpy.typing as npt
import zarr
from zarr.core.group import Group, GroupMetadata

from atomea.data import OptionalSliceSpec
from atomea.stores import DiskFormat
from atomea.stores.arrays import ArrayStore


class ZarrArrayStore(ArrayStore):
    """
    Zarr-based array store where each logical path maps to a Zarr array.

    Paths use '/' to denote nested logical structure.
    """

    def __init__(
        self,
        path: str,
        disk_format: DiskFormat = DiskFormat.ZARR,
        mode: str = "r",
        **kwargs: Any,
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
        self._store: Group = zarr.open_group(store=path, mode=mode, **kwargs)
        super().__init__(path, disk_format=disk_format, **kwargs)

    def create(
        self,
        path: str,
        shape: tuple[int, ...],
        *args: Any,
        overwrite: bool = False,
        dtype: npt.DTypeLike | None = None,
        chunks: tuple[int, ...] | Literal["auto"] = "auto",
        **kwargs: Any,
    ) -> Any:
        """
        Pre-allocate an array with the given shape and dtype.

        Args:
            path: hierarchical key, e.g. 'coords'.
            shape: full shape of the array.
            overwrite: if True, delete existing before create; otherwise error if exists.
            dtype: numpy-compatible dtype.
            chunks: chunk shape or 'auto'.
        """
        group_path, _ = path.rsplit("/", 1) if "/" in path else ("", path)
        self._store.create_hierarchy({group_path: GroupMetadata()})
        z = self._store.create_array(
            name=path,
            shape=shape,
            dtype=dtype,
            chunks=chunks,
            overwrite=overwrite,
            **kwargs,
        )
        return z

    def write(
        self,
        path: str,
        data: npt.NDArray[np.generic],
        view: OptionalSliceSpec = None,
        dtype: npt.DTypeLike | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Write data to an array, whole or sliced. Will create the array if it does not
        exist.

        Args:
            path: hierarchical key of the array.
            array: data to write.
            view: if None, overwrite entire array (requires create or existing),
                    else write into the specified slice region.
            overwrite: if True and slices is None, replaces existing array definition.
        """
        z = self._store.get(path=path)
        if not dtype:
            dtype = data.dtype
        if z is None:
            z = self.create(path, data.shape, dtype=dtype, **kwargs)
            z[:] = data
        else:
            z.set_basic_selection(view, data, **kwargs)  # type: ignore

    def append(
        self, path: str, data: npt.NDArray[np.generic], *args: Any, **kwargs: Any
    ) -> None:
        """
        Append data along the first axis to an existing Zarr array;
        creates the array with an unlimited first dimension if it does not exist.
        """
        arr = self.read(path)
        arr.append(data, **kwargs)  # type: ignore

    def read(
        self,
        path: str,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> npt.NDArray[np.generic] | None:
        """
        Read the array from Zarr, optionally returning a subset efficiently.
        Returns None if the path does not exist.
        """
        z = self._store.get(path)
        if view is None:
            return z[:]  # type: ignore
        if isinstance(view, dict):
            idx = [view(None)] * z.ndim  # type: ignore
            for axis, sl in view.items():
                idx[axis] = sl
            return z[tuple(idx)]  # type: ignore
        return z[view]  # type: ignore

    def available(self) -> list[str]:
        """
        List all stored Zarr array paths (joined by '/').
        """
        keys = list(self._store.array_keys())
        return keys
