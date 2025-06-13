"""Type system for the atomea Data descriptor.

This module provides a comprehensive type system for scientific data storage and access,
supporting both NumPy arrays of various precisions and Polars DataFrames. The type system
enables type-safe data operations while maintaining flexibility for different data types
commonly used in atomistic simulations and analysis.

The module defines:
    - Generic type variables for flexible data handling
    - Slice specifications for multidimensional array access
    - NumPy array type aliases for different numeric precisions
    - Optional versions of all data types for missing data handling
    - DataFrame types for tabular data

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

See Also:
    atomea.data.Data: The main data descriptor class
    atomea.data.accessors: Array and table accessor classes
    atomea.stores: Storage backend implementations
"""

from ._generic import T, SliceSpec, ValueOrSlice
from ._array import (
    Float64,
    Float32,
    Int64,
    Int32,
    Int8,
    UInt64,
    UInt32,
    UInt8,
    Str,
    OptionalFloat64,
    OptionalFloat32,
    OptionalInt64,
    OptionalInt32,
    OptionalUInt8,
    OptionalStr,
    DTYPE_TO_NUMPY,
)
from ._dataframe import DataFrame, OptionalDataFrame

__all__ = [
    "T",
    "SliceSpec",
    "ValueOrSlice",
    "Float64",
    "Float32",
    "Int64",
    "Int32",
    "Int8",
    "UInt64",
    "UInt32",
    "UInt8",
    "Str",
    "OptionalFloat64",
    "OptionalFloat32",
    "OptionalInt64",
    "OptionalInt32",
    "OptionalUInt8",
    "OptionalStr",
    "DataFrame",
    "OptionalDataFrame",
    "DTYPE_TO_NUMPY",
]
