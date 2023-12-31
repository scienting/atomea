from typing import Any

from collections.abc import Iterable

import numpy.typing as npt
import zarr

DTYPE_MAP = {
    "int8": "i1",
    "int16": "i2",
    "int32": "i4",
    "int64": "i8",
    "uint8": "u1",
    "uint16": "u2",
    "uint32": "u4",
    "uint64": "u8",
    "float32": "f4",
    "float64": "f8",
    "utf-8": "U64",
}


def initialize(
    path: str, shape: tuple[int, ...], dtype: str, **kwargs: dict[str, Any]
) -> zarr.core.Array:
    """Thin wrapper around `zarr.creation.create`."""
    if not path.endswith(".zarr"):
        path += ".zarr"
    return zarr.creation.create(shape, dtype=DTYPE_MAP[dtype], path=path, **kwargs)


def write(
    arr: zarr.core.Array,
    data: npt.NDArray[Any],
    indices: int | Iterable[int],
    *args: Any,
    **kwargs: dict[str, Any],
) -> zarr.core.Array:
    arr[indices] = data
    return arr
