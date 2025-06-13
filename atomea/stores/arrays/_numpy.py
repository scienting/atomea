import os

import numpy as np
import numpy.typing as npt

from atomea.stores import DiskFormat
from atomea.stores.arrays import ArrayStore


class NumpyArrayStore(ArrayStore):
    """
    In-memory array store using NumPy arrays in a flat dict keyed by hierarchical paths.

    Paths use '/' to denote nested logical structure but are stored as flat keys.
    """

    def __init__(self, disk_format: DiskFormat = DiskFormat.NPZ) -> None:
        self._store: dict[str, npt.NDArray[np.generic]] = {}
        super().__init__(disk_format)

    def write(self, path: str, array: npt.NDArray[np.generic]) -> None:
        """
        Write or overwrite the array at the given path.
        """
        self._store[path] = np.array(array, copy=True)

    def append(self, path: str, array: npt.NDArray[np.generic]) -> None:
        """
        Append data along the first axis to an existing array at path;
        if no array exists, creates one.
        """
        arr = np.array(array, copy=False)
        if path in self._store:
            existing = self._store[path]
            try:
                self._store[path] = np.concatenate((existing, arr), axis=0)
            except ValueError as e:
                raise ValueError(
                    f"Cannot append array with shape {arr.shape} to existing "
                    f"array of shape {existing.shape}"
                ) from e
        else:
            self.write(path, arr)

    def read(
        self,
        path: str,
        slices: tuple[slice, ...] | dict[int, slice | tuple[slice, ...]] | None = None,
    ) -> npt.NDArray[np.generic] | None:
        """
        Read the array at path, optionally returning a subset.

        Args:
            path: the key of the stored array.
            slices: either a tuple of slice objects for each axis,
                or a dict mapping axis index to a tuple of slice objects.
                If None, returns the full array.

        Returns:
            The requested ndarray, or None if path not found.
        """
        data = self._store.get(path)
        if data is None:
            return None
        if slices is None:
            return data
        if isinstance(slices, dict):  # axis-wise slicing
            full = [slice(None)] * data.ndim
            for ax, sl in slices.items():
                full[ax] = sl
            return data[tuple(full)]
        # tuple of slices
        return data[slices]

    def available(self) -> list[str]:
        """
        List all stored paths.
        """
        return list(self._store.keys())

    def _dump_npy(self, prefix: str = "") -> None:
        for key, arr in self._store:
            path = os.path.join(prefix, key + ".npy")
            np.save(path, arr)

    def _dump_npz(self, prefix: str = "") -> None:
        if not prefix:
            prefix = "arrays"
        fn = prefix if prefix.endswith(".npz") else prefix + ".npz"

        parent = os.path.dirname(fn)
        if parent and not os.path.isdir(parent):
            os.makedirs(parent, exist_ok=True)

        # build a dict with “safe” names
        savez_dict = {key.replace("/", "_"): arr for key, arr in self._store.items()}
        np.savez(fn, *savez_dict)

    def dump(self, prefix: str = "") -> None:
        if self.disk_format == DiskFormat.NPY:
            self._dump_npy(prefix)
        elif self.disk_format == DiskFormat.NPZ:
            self._dump_npz(prefix)
        else:
            raise ValueError("DiskFormat of {} not supported", self.disk_format)
