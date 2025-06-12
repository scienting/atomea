import numpy as np
import numpy.typing as npt

from atomea.stores.arrays import ArrayStore


class NumpyArrayStore(ArrayStore):
    """
    In-memory array store using NumPy arrays in a flat dict keyed by hierarchical paths.

    Paths use '/' to denote nested logical structure but are stored as flat keys.
    """

    def __init__(self) -> None:
        self._store: dict[str, npt.NDArray[np.generic]] = {}

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
        slices: tuple[slice, ...] | dict[int, tuple[slice, ...]] | None = None,
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
