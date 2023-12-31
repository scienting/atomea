from typing import Any

from collections.abc import Iterable

import pyarrow as pa
import pyarrow.parquet as pq

from .utils import apply_dtype_map

dtype_map = {
    "int8": "int8",
    "int16": "int16",
    "int32": "int32",
    "int64": "int64",
    "uint8": "uint8",
    "uint16": "uint16",
    "uint32": "uint32",
    "uint64": "uint64",
    "float32": "float32",
    "float64": "float64",
    "utf-8": "string",
}


def write(
    path: str,
    data: dict[str, Any],
    schema_fields: Iterable[Iterable[Any]],
    schema_metadata: dict[str, Any] | None = None,
    **kwargs: dict[str, Any],
) -> pa.Table:
    """Initialize `pyarrow.table` and save parquet file.

    Args:
        path: Path to parquet file to save.
        data: Columns (keys) and values to store in a parquet file.
        schema_fields: TODO: Document this.
        **kwargs: Passed into `pa.table` and `pa.write_table`.
    """
    if not path.endswith(".parquet"):
        path += ".parquet"
    schema_fields = apply_dtype_map(schema_fields, dtype_map)
    if schema_metadata is None:
        schema_metadata = {}
    schema = pa.schema(schema_fields, metadata=schema_metadata)
    arrays = [pa.array(data[k], type=schema.field_by_name(k).type) for k in data]
    table = pa.table(arrays, schema=schema, **kwargs)
    pq.write_table(table=table, where=path, **kwargs)
    return table
