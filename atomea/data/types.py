from typing import TypeAlias, TypeVar

T = TypeVar("T")
"""Generic type variable for Data descriptor parametrization.

This creates a generic type variable that can represent any supported data type.
It enables the `Data[T]` descriptor to be parameterized with specific types,
providing compile-time type safety and runtime type consistency.

When used in a `Data[T]` declaration, `T` gets bound to the specific type alias,
allowing type checkers to understand the expected data type and provide
appropriate code completion and error detection.

Example:
    Basic usage::

        # T gets bound to Float64 type alias
        field: Data[Float64] = Data[Float64](dtype=np.float64, ...)

        # Type checker knows this returns Float64 | None
        data = obj.field.value

    Generic function usage::

        def process_data(field: Data[T]) -> T | None:
            return field.value  # Type checker knows return type

Note:
    The type variable is unbounded, meaning it can represent any type.
    However, in practice it should only be used with the type aliases
    defined in this module for consistent behavior.
"""

Slice = tuple[slice | int | list[int] | tuple[int, ...], ...]

SliceSpec = Slice | dict[int, Slice]
"""Type specification for multidimensional array slicing operations.

Defines the allowed formats for specifying how to slice multidimensional arrays
in both read and write operations. Supports both simple tuple-based slicing
and advanced axis-specific slicing for complex data access patterns.

The type union supports two slicing approaches:

1. **Tuple format** (`tuple[slice | int | list[int] | tuple[int, ...], ...]`):
   - Direct specification of slice objects for each dimension
   - Maps directly to standard NumPy indexing syntax
   - Most common and intuitive approach

2. **Dictionary format** (`dict[int, tuple[slice, ...]]`):
   - Maps axis indices to slice specifications
   - Allows sparse slicing (only specify axes you want to slice)
   - Useful for high-dimensional arrays where you only want to slice specific axes

Example:
    Tuple format slicing::

        # Slice first 10 elements along axis 0, all elements along axis 1
        slice_spec = (slice(0, 10), slice(None))
        data = interface(slices=view_spec)

        # Equivalent to: array[0:10, :]

    Dictionary format slicing::

        # Slice axis 0 and axis 2, leave axis 1 unchanged
        slice_spec = {
            0: (slice(0, 10),),  # First 10 elements on axis 0
            2: (slice(None, None, 2),)  # Every 2nd element on axis 2
        }
        data = interface(slices=view_spec)

        # For 3D array, equivalent to: array[0:10, :, ::2]

    Direct interface usage::

        # Using bracket notation (converted to SliceSpec internally)
        first_10_atoms = coordinates[0:10, :]
        x_coordinates = coordinates[:, 0]
        every_other_timestep = coordinates[:, :, ::2]

Note:
    The slice specifications are used internally by the interface classes
    and are typically not created directly by users. The bracket notation
    on interfaces (`interface[0:10]`) automatically converts to appropriate
    SliceSpec formats.
"""

OptionalSliceSpec = SliceSpec | None


ValueOrSlice: TypeAlias = T | tuple[T, SliceSpec]
"""Type specification for data values that support partial writing.

Represents data that can be either a complete value or a value with associated
slice specification for partial array updates. This type is primarily used
in setter operations where you might want to update only a portion of an
existing array rather than replacing the entire array.

The type union supports two writing modes:

1. **Complete replacement** (`T`):
   - Replaces the entire stored array/data
   - Most common usage pattern
   - Simpler and more straightforward

2. **Partial update** (`tuple[T, SliceSpec]`):
   - Updates only the specified slice region
   - Preserves existing data outside the slice
   - Efficient for large arrays with small updates

Type Parameters:
    T: The underlying data type (e.g., Float64, DataFrame)

Example:
    Complete array replacement::

        # Replace entire coordinate array
        new_coords = np.random.rand(100, 3)
        ensemble.coordinates = new_coords

    Partial array updates::

        # Update only first 10 atoms' coordinates
        partial_coords = np.random.rand(10, 3)
        slice_spec = (slice(0, 10), slice(None))
        ensemble.coordinates = (partial_coords, slice_spec)

        # Update only x-coordinates (axis 1, index 0)
        x_coords = np.random.rand(100)
        ensemble.coordinates = (x_coords.reshape(-1, 1), (slice(None), slice(0, 1)))

    Usage in Data descriptor::

        def __set__(self, obj, value: ValueOrSlice[T]) -> None:
            if isinstance(value, tuple) and len(value) == 2:
                data, slices = value
                # Partial update
                store.write(path, data, slices=views)
            else:
                # Complete replacement
                store.write(path, value, slices=None)

Warning:
    Partial updates require that the target array already exists and has
    compatible dimensions. The slice specification must be valid for the
    existing array shape, or a runtime error will occur.

See Also:
    SliceSpec: Slice specification format
    atomea.data.Data.__set__: Implementation of partial writing
    atomea.stores: Storage backend partial write support
"""
