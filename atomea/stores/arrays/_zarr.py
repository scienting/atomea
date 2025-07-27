# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any, Literal

from collections.abc import Iterator
from pathlib import Path

import numpy as np
import numpy.typing as npt
import zarr
from zarr.core.group import Group, GroupMetadata

from atomea.data import OptionalSliceSpec
from atomea.helpers import chunker
from atomea.stores import DiskFormat
from atomea.stores.arrays import ArrayStore


class ZarrArrayStore(ArrayStore):
    """
    Zarr-based array store where each logical path maps to a Zarr array.

    Paths use '/' to denote nested logical structure.
    """

    def __init__(
        self,
        path: Path | str,
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.ZARR,
        **kwargs: Any,
    ) -> None:
        assert disk_format == DiskFormat.ZARR
        self._store: Group = zarr.open_group(store=str(path), mode=mode, **kwargs)
        super().__init__(path, mode=mode, disk_format=disk_format, **kwargs)

    def create(
        self,
        path: Path | str,
        shape: tuple[int, ...],
        overwrite: bool = False,
        dtype: npt.DTypeLike | None = None,
        chunks: tuple[int, ...] | Literal["auto"] = "auto",
        **kwargs: Any,
    ) -> Any:
        if self.mode == "r":
            raise ValueError("Cannot create when in 'r' mode")
        path = str(path)
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
        path: Path | str,
        data: npt.NDArray[np.generic],
        view: OptionalSliceSpec = None,
        dtype: npt.DTypeLike | None = None,
        **kwargs: Any,
    ) -> None:
        if self.mode == "r":
            raise ValueError("Cannot write when in 'r' mode")
        z = self._store.get(path=str(path))
        if not dtype:
            dtype = data.dtype
        if z is None:
            if self.mode in ("r+", "w-"):
                raise ValueError(f"Cannot create data when in '{self.mode}' mode")
            z = self.create(path, data.shape, dtype=dtype, **kwargs)
            z[:] = data
        else:
            z.set_basic_selection(view, data, **kwargs)  # type: ignore

    def append(
        self, path: Path | str, data: npt.NDArray[np.generic], *args: Any, **kwargs: Any
    ) -> None:
        if self.mode in ("r", "w-"):
            raise ValueError(f"Cannot append when in '{self.mode}' mode")
        arr = self.read(str(path))
        arr.append(data, **kwargs)  # type: ignore

    def get(
        self,
        path: Path | str,
        **kwargs: Any,
    ) -> zarr.Array | None:
        z = self._store.get(str(path))
        assert not isinstance(z, zarr.Group)
        return z

    def read(
        self,
        path: Path | str,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> npt.NDArray[np.generic] | None:
        z = self.get(path)
        if view is None:
            return z[:]  # type: ignore
        if isinstance(view, dict):
            idx = [view(None)] * z.ndim  # type: ignore
            for axis, sl in view.items():
                idx[axis] = sl
            return z[tuple(idx)]  # type: ignore
        return z[view]  # type: ignore

    def available(self) -> list[str]:
        keys = list(self._store.array_keys())
        return keys

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
            yield z.get_orthogonal_selection(_view)
