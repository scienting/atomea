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
