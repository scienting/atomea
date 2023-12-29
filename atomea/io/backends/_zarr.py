from typing import Any

import numpy.typing as npt
import zarr

dtype_map = {
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
}


def array_initialize(path: str, shape: tuple[int, ...], dtype: str, **kwargs):
    """Thin wrapper around `zarr.creation.create`."""
    if not path.endswith(".zarr"):
        path += ".zarr"
    return zarr.creation.create(shape, dtype=dtype_map[dtype], path=path, **kwargs)


def array_write(
    path: str,
    arr: npt.NDArray[Any],
    mode: str = "a",
    dtype: str = "float64",
    **kwargs: dict[str, Any],
) -> zarr.core.Array:
    """Thin wrapper around `zarr.open`."""
    if not path.endswith(".zarr"):
        path += ".zarr"
    dtype = dtype_map[dtype]
    z = zarr.open(path, mode=mode, shape=arr.shape, dtype=dtype, **kwargs)
    z[:] = arr
    return z
