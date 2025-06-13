from enum import Enum


class StoreKind(Enum):
    ARRAY = 0
    TABLE = 1


class DiskFormat(Enum):
    NONE = 0

    # Arrays
    ZARR = 10
    NPY = 11
    NPZ = 12

    # TABLES
    CSV = 20
    IPC = 21
    EXCEL = 22
    JSON = 23
    PARQUET = 24
    AVRO = 25
    DELTA = 26
    ICEBERG = 27


ArrayDiskFormats = (DiskFormat.ZARR, DiskFormat.NPY, DiskFormat.NPZ)
