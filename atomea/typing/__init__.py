"""Type system for the atomea Data descriptor.

This module provides a comprehensive type system for scientific data storage and access,
supporting both NumPy arrays of various precisions and Polars DataFrames. The type system
enables type-safe data operations while maintaining flexibility for different data types
commonly used in atomistic simulations and analysis.

The module defines:

- NumPy array type aliases for different numeric precisions;
- Optional versions of all data types for missing data handling;
- DataFrame types for tabular data.

Example:
    Creating type-safe data fields::

        # 64-bit float coordinates
        coordinates: Data[Float64] = Data[Float64](dtype=np.float64, meta=Metadata(...))

        # 8-bit unsigned integers for atomic numbers
        atom_numbers: Data[Uint8] = Data[Uint8](dtype=np.uint8, meta=Metadata(...))

        # DataFrame for tabular results
        results: Data[DataFrame] = Data[DataFrame](dtype=None, meta=Metadata(...))

    Type-safe access::

        # Type checkers understand these return the correct types
        coords: Float64 | None = ensemble.coordinates.value
        atoms: Uint8 | None = ensemble.atom_numbers.value
        data: pl.DataFrame | None = ensemble.results.value
"""

from ._array import (
    Float64,
    OptionalFloat64,
    Float32,
    OptionalFloat32,
    Int64,
    OptionalInt64,
    Int32,
    OptionalInt32,
    Int16,
    OptionalInt16,
    Int8,
    OptionalInt8,
    UInt64,
    OptionalUInt64,
    UInt32,
    OptionalUInt32,
    UInt16,
    OptionalUInt16,
    UInt8,
    OptionalUInt8,
    Str,
    OptionalStr,
)
from ._dataframe import DataFrame, OptionalDataFrame
from .passable import PassableData, OptionalPassableData

__all__ = [
    "Float64",
    "OptionalFloat64",
    "Float32",
    "OptionalFloat32",
    "Int64",
    "OptionalInt64",
    "Int32",
    "OptionalInt32",
    "Int16",
    "OptionalInt16",
    "Int8",
    "OptionalInt8",
    "UInt64",
    "OptionalUInt64",
    "UInt32",
    "OptionalUInt32",
    "UInt16",
    "OptionalUInt16",
    "UInt8",
    "OptionalUInt8",
    "Str",
    "OptionalStr",
    "DataFrame",
    "OptionalDataFrame",
    "PassableData",
    "OptionalPassableData",
]
