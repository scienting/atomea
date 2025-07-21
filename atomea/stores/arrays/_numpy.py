# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

import os
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec
from atomea.helpers import chunker
from atomea.stores import DiskFormat
from atomea.stores.arrays import ArrayStore


class NumpyArrayStore(ArrayStore):
    """
    In-memory array store using NumPy arrays in a flat dict keyed by hierarchical paths.

    Paths use '/' to denote nested logical structure but are stored as flat keys.
    """

    def __init__(
        self,
        path: Path | str,
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.NPZ,
        **kwargs: Any,
    ) -> None:
        self._store: dict[str, npt.NDArray[np.generic]] = {}
        super().__init__(path, mode=mode, disk_format=disk_format, **kwargs)

    def create(
        self,
        path: Path | str,
        shape: tuple[int, ...],
        overwrite: bool = False,
        dtype: npt.DTypeLike = np.dtype(np.float64),
        fill: np.generic | None = None,
        **kwargs: Any,
    ) -> Any:
        if self.mode == "r":
            raise ValueError("Cannot write when in 'r' mode")
        path = str(path)
        if path in self._store and not overwrite:
            raise RuntimeError(f"{path} already exists and overwrite is False!")
        if fill:
            self._store[path] = np.full(shape, fill, dtype=dtype)
        else:
            self._store[path] = np.empty(shape, dtype=dtype)

    def write(
        self,
        path: Path | str,
        data: npt.NDArray[np.generic],
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        if self.mode == "r":
            raise ValueError("Cannot write when in 'r' mode")
        path = str(path)
        if path in self._store.keys():
            if self.mode in ("r+", "w-"):
                raise ValueError(f"Cannot create data when in '{self.mode}' mode")
        self._store[path] = np.array(data, copy=True)

    def append(
        self, path: Path | str, data: npt.NDArray[np.generic], *args: Any, **kwargs: Any
    ) -> None:
        if self.mode in ("r", "w-"):
            raise ValueError(f"Cannot append when in '{self.mode}' mode")
        path = str(path)
        arr = np.array(data, copy=False)
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

    def get(
        self,
        path: Path | str,
        **kwargs: Any,
    ) -> npt.NDArray[np.generic] | None:
        path = str(path)
        data = self._store.get(path)
        return data

    def read(
        self, path: Path | str, view: OptionalSliceSpec = None, **kwargs: Any
    ) -> npt.NDArray[np.generic] | None:
        data = self.get(path)
        if data is None:
            return None
        if view is None:
            return data
        if isinstance(view, dict):  # axis-wise slicing
            full = [slice(None)] * data.ndim
            for ax, sl in view.items():
                full[ax] = sl
            return data[tuple(full)]
        # tuple of view
        return data[view]

    def available(self) -> list[str]:
        return list(self._store.keys())

    def _flush_npy(self, prefix: str = "") -> None:
        for key, arr in self._store:
            path = os.path.join(prefix, key + ".npy")
            np.save(path, arr)

    def _flush_npz(self, prefix: str = "") -> None:
        if not prefix:
            prefix = "arrays"
        fn = prefix if prefix.endswith(".npz") else prefix + ".npz"

        parent = os.path.dirname(fn)
        if parent and not os.path.isdir(parent):
            os.makedirs(parent, exist_ok=True)

        # build a dict with “safe” names
        savez_dict = {key.replace("/", "_"): arr for key, arr in self._store.items()}
        np.savez(fn, *savez_dict)

    def flush(self, **kwargs: Any) -> None:
        path = str(self.path)
        if self.disk_format == DiskFormat.NPY:
            self._flush_npy(path)
        elif self.disk_format == DiskFormat.NPZ:
            self._flush_npz(path)
        else:
            raise ValueError(f"DiskFormat of {self.disk_format} not supported")

    def iter(
        self,
        path: Path | str,
        elements: OptionalSliceSpec = None,
        view: OptionalSliceSpec = None,
        chunk_size: int = 1,
        **kwargs: Any,
    ) -> Iterator[Any]:
        """Yield chunks of data instead of reading all into memory.

        Args:
            path:
            view: View for all but the first dimension.
            chunk: Number of data points of the first axis to yield.
        """
        z = self.get(path)
        if z is None:
            return None
        n_items = z.shape[0]
        for chunk in chunker(chunk_size, n_items, elements):
            if view is None:
                _view = (chunk,)
            else:
                _view = (chunk, *view)  # type: ignore
            yield z[_view]
