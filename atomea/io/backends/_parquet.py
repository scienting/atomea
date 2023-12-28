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
    "utf8": "string",
}


def tabular_parquet(
    path: str,
    data: dict[str, Any],
    schema_fields: Iterable[Iterable[Any]],
    **kwargs: dict[str, Any],
) -> pa.Table:
    """Initialize `pyarrow.table` and save parquet file."""
    if not path.endswith(".parquet"):
        path += ".parquet"
    schema_fields = apply_dtype_map(schema_fields, dtype_map)
    schema = pa.schema(schema_fields)

    arrays = [pa.chunked_array(v) for k, v in data.items()]
    table = pa.table(arrays, schema=schema, **kwargs)

    pq.write_table(table=table, where=path, **kwargs)
    return table
