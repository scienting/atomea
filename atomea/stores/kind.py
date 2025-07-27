# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from enum import Enum, auto


class StoreKind(Enum):
    """The kind of storage backend."""

    ARRAY = auto()
    TABLE = auto()


class DiskFormat(Enum):
    """Format of data when flushed to disk."""

    NONE = auto()

    # Arrays
    ZARR = auto()
    NPY = auto()
    NPZ = auto()

    # TABLES
    CSV = auto()
    IPC = auto()
    EXCEL = auto()
    JSON = auto()
    PARQUET = auto()
    AVRO = auto()
    DELTA = auto()
    ICEBERG = auto()


ArrayDiskFormats = (DiskFormat.ZARR, DiskFormat.NPY, DiskFormat.NPZ)
