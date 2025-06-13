from typing import TypeAlias

import numpy as np
import numpy.typing as npt

Float64: TypeAlias = npt.NDArray[np.float64]
"""64-bit double-precision floating-point array type.

Represents NumPy arrays containing 64-bit floating-point numbers (double precision).
This is the default floating-point type in NumPy and provides the highest precision
for floating-point calculations, making it suitable for scientific computations
where numerical accuracy is critical.

Characteristics:
    - Precision: ~15-17 decimal digits
    - Range: ±1.8 × 10^308
    - Memory: 8 bytes per element
    - IEEE 754 double precision standard

Typical Use Cases:
    - Atomic coordinates (x, y, z positions)
    - Energy values (electronic, kinetic, potential)
    - Physical properties requiring high precision
    - Intermediate calculations in scientific algorithms
    - Reference data and benchmarks

Example:
    Field definition::

        coordinates: Data[Float64] = Data[Float64](
            dtype=np.float64,
            meta=Metadata(
                description="Atomic coordinates in Angstroms",
                store=StoreKind.ARRAY
            )
        )

    Creating data::

        # 100 atoms with x,y,z coordinates
        coords = np.random.rand(100, 3).astype(np.float64)
        ensemble.coordinates = coords

    Type checking::

        # Type checker knows this is npt.NDArray[np.float64]
        coord_data: Float64 = ensemble.coordinates.value

Note:
    This type alias does not include None. For optional arrays, use OptionalFloat64.
    The high precision comes with increased memory usage compared to Float32.

See Also:
    Float32: Single-precision alternative for memory efficiency
    OptionalFloat64: Optional version including None
    numpy.float64: NumPy documentation for float64 type
"""

Float32: TypeAlias = npt.NDArray[np.float32]
"""32-bit single-precision floating-point array type.

Represents NumPy arrays containing 32-bit floating-point numbers (single precision).
Provides a good balance between precision and memory efficiency, suitable for
applications where memory usage is a concern but reasonable floating-point
precision is still required.

Characteristics:
    - Precision: ~6-9 decimal digits
    - Range: ±3.4 × 10^38
    - Memory: 4 bytes per element (half of Float64)
    - IEEE 754 single precision standard

Typical Use Cases:
    - Large datasets where memory is constrained
    - Intermediate calculations with acceptable precision loss
    - Graphics and visualization data
    - Machine learning features
    - Approximate physical quantities

Example:
    Memory-efficient coordinates::

        coordinates: Data[Float32] = Data[Float32](
            dtype=np.float32,
            meta=Metadata(
                description="Coordinates (memory optimized)",
                store=StoreKind.ARRAY
            )
        )

    Large dataset storage::

        # 1 million atoms - saves 12MB compared to Float64
        large_coords = np.random.rand(1000000, 3).astype(np.float32)
        ensemble.coordinates = large_coords

Performance Considerations:
    - Uses 50% less memory than Float64
    - May be faster on some hardware due to cache efficiency
    - Precision loss may accumulate in iterative calculations
    - Consider precision requirements vs. memory constraints

Warning:
    Single precision may not be sufficient for high-accuracy scientific
    calculations. Evaluate precision requirements before choosing Float32
    over Float64 for critical computations.

See Also:
    Float64: Double-precision alternative for higher accuracy
    OptionalFloat32: Optional version including None
    numpy.float32: NumPy documentation for float32 type
"""


Int64: TypeAlias = npt.NDArray[np.int64]
"""64-bit signed integer array type.

Represents NumPy arrays containing 64-bit signed integers (long integers).
Provides the largest range for integer values in standard NumPy types,
suitable for applications requiring very large integer values or when
maximum precision is needed for integer calculations.

Characteristics:
    - Range: -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807 (±9.2 × 10^18)
    - Memory: 8 bytes per element
    - Signed: Can represent negative values
    - Default integer type on 64-bit systems

Typical Use Cases:
    - Very large molecular systems (billions of atoms)
    - High-precision time stamps or step counters
    - File offsets and memory addresses
    - Large-scale simulation indices
    - Unique identifiers requiring large ranges

Example:
    Large system indices::

        atom_ids: Data[Int64] = Data[Int64](
            dtype=np.int64,
            meta=Metadata(
                description="Global atom identifiers",
                store=StoreKind.ARRAY
            )
        )

    Timestamp data::

        # Nanosecond precision timestamps
        timestamps = np.array([1234567890123456789, 1234567890123456790], dtype=np.int64)
        simulation.timestamps = timestamps

    Large index ranges::

        # For systems with billions of elements
        global_indices = np.arange(0, 2**40, dtype=np.int64)

Memory Considerations:
    - Uses twice the memory of Int32
    - Consider whether the full range is actually needed
    - May impact cache performance due to larger size

Performance Notes:
    - Native integer size on 64-bit systems
    - No overflow concerns for most practical applications
    - Slower than smaller integer types on memory-constrained systems

See Also:
    Int32: More memory-efficient 32-bit version
    OptionalInt64: Optional version including None
    numpy.int64: NumPy documentation for int64 type
"""


Int32: TypeAlias = npt.NDArray[np.int32]
"""32-bit signed integer array type.

Represents NumPy arrays containing 32-bit signed integers. Provides a good
balance between range and memory efficiency for integer data that doesn't
require the full range of 64-bit integers.

Characteristics:
    - Range: -2,147,483,648 to 2,147,483,647 (approximately ±2.1 billion)
    - Memory: 4 bytes per element
    - Signed: Can represent negative values
    - No fractional component

Typical Use Cases:
    - Atom indices and connectivity information
    - Bond lists and topology data
    - Timestep indices in trajectories
    - Grid indices for spatial decomposition
    - Reference indices and mappings

Example:
    Bond connectivity::

        bonds: Data[Int32] = Data[Int32](
            dtype=np.int32,
            meta=Metadata(
                description="Bond connectivity pairs",
                store=StoreKind.ARRAY
            )
        )

    Creating bond data::

        # Bonds between atoms (atom1_idx, atom2_idx)
        bond_pairs = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int32)
        topology.bonds = bond_pairs

    Index arrays::

        # Indices of selected atoms
        selected_atoms = np.array([0, 5, 10, 15, 20], dtype=np.int32)

Range Considerations:
    - Sufficient for most molecular systems (millions of atoms)
    - Consider Int64 for very large systems or when indexing huge arrays
    - Negative values allowed (useful for special markers/flags)

Performance Notes:
    - More memory efficient than Int64
    - Standard integer size on many 32-bit systems
    - Good balance of range and efficiency

See Also:
    Int64: 64-bit version for larger ranges
    Uint8: Unsigned 8-bit for small positive integers
    OptionalInt32: Optional version including None
    numpy.int32: NumPy documentation for int32 type
"""

Int8: TypeAlias = npt.NDArray[np.int8]

UInt64: TypeAlias = npt.NDArray[np.uint64]

UInt32: TypeAlias = npt.NDArray[np.uint32]

UInt8: TypeAlias = npt.NDArray[np.uint8]
"""8-bit unsigned integer array type.

Represents NumPy arrays containing 8-bit unsigned integers. This is the most
memory-efficient integer type, suitable for small positive integer values
such as atomic numbers, small counts, or categorical data.

Characteristics:
    - Range: 0 to 255 (256 unique values)
    - Memory: 1 byte per element (most efficient)
    - Unsigned: Only non-negative values
    - Commonly used for compact data representation

Typical Use Cases:
    - Atomic numbers (1-118 for all elements)
    - Small categorical identifiers
    - Boolean-like flags (0/1)
    - Color values and image data
    - Compact enumeration indices

Example:
    Atomic numbers::

        atom_z: Data[Uint8] = Data[Uint8](
            dtype=np.uint8,
            meta=Metadata(
                description="Atomic numbers",
                store=StoreKind.ARRAY
            )
        )

    Creating atomic number data::

        # H, H, O for water molecule
        atomic_numbers = np.array([1, 1, 8], dtype=np.uint8)
        molecule.atom_z = atomic_numbers

    Element mapping::

        # All periodic table elements fit in uint8
        all_elements = np.arange(1, 119, dtype=np.uint8)  # H to Og

Memory Efficiency:
    - 8x more memory efficient than Int64
    - 4x more memory efficient than Int32
    - Excellent for large arrays of small values
    - Significant memory savings for molecular systems

Range Limitations:
    - Cannot represent negative values
    - Maximum value is 255
    - Consider larger integer types if range is insufficient

Use Cases in Chemistry:
    - Perfect for atomic numbers (max 118)
    - Residue types in proteins
    - Bond orders (1, 2, 3)
    - Small multiplicity values

See Also:
    Int32: Signed integers with larger range
    OptionalUint8: Optional version including None
    numpy.uint8: NumPy documentation for uint8 type
"""

Str: TypeAlias = npt.NDArray[np.str_]
"""NumPy string array type for text data.

Represents NumPy arrays containing string values. These arrays store
fixed-length Unicode strings and are useful for categorical text data
such as element symbols, residue names, or other textual identifiers
in scientific datasets.

Characteristics:
    - Fixed-length Unicode strings
    - Memory efficient for repeated short strings
    - NumPy array operations supported
    - Vectorized string operations available

Typical Use Cases:
    - Chemical element symbols ('H', 'He', 'Li', etc.)
    - Residue names in proteins ('ALA', 'GLY', 'SER', etc.)
    - Atom types in force fields ('CT', 'HC', 'OH', etc.)
    - Categorical labels and identifiers
    - File names or path components

Example:
    Element symbols::

        atom_symbols: Data[Str] = Data[Str](
            dtype=np.str_,
            meta=Metadata(
                description="Chemical element symbols",
                store=StoreKind.ARRAY
            )
        )

    Creating symbol data::

        # Water molecule symbols
        symbols = np.array(['H', 'H', 'O'], dtype=np.str_)
        molecule.atom_symbols = symbols

    Protein residues::

        # Protein sequence
        residues = np.array(['MET', 'ALA', 'GLY', 'SER'], dtype=np.str_)
        protein.residue_names = residues

String Operations:
    - Vectorized comparison: `symbols == 'H'`
    - String methods: `np.char.upper(symbols)`
    - Pattern matching and searching
    - Concatenation and manipulation

Memory Considerations:
    - All strings in array have same maximum length
    - Shorter strings are padded with null characters
    - Consider Python object arrays for variable-length strings
    - More memory efficient than object arrays for short, similar strings

Performance Notes:
    - Fast for fixed-length string operations
    - NumPy vectorization benefits
    - Suitable for categorical data analysis
    - Less flexible than Python string objects

See Also:
    OptionalStr: Optional version including None
    numpy.str_: NumPy documentation for string arrays
    numpy.char: NumPy string operation functions
"""


OptionalFloat64: TypeAlias = Float64 | None
"""Optional 64-bit floating-point array type.

Extends Float64 to include None, representing arrays that may be missing
or not yet initialized. This is the most commonly used type for optional
floating-point data in scientific applications.

Type Definition:
    `npt.NDArray[np.float64] | None`

Typical Use Cases:
    - Data fields that may not be available for all systems
    - Optional computed properties
    - Conditional analysis results
    - Lazy-loaded expensive calculations

Example:
    Optional energy field::

        energy: Data[OptionalFloat64] = Data[OptionalFloat64](
            dtype=np.float64,
            meta=Metadata(
                description="Total energy (computed on demand)",
                store=StoreKind.ARRAY
            ),
            default=None
        )

    Handling optional data::

        energy_data = ensemble.energy.value
        if energy_data is not None:
            max_energy = np.max(energy_data)
        else:
            print("Energy not computed yet")

    Type-safe processing::

        def process_energy(data: OptionalFloat64) -> float | None:
            return np.mean(data) if data is not None else None

See Also:
    Float64: Non-optional version
    Data.default: Setting default values for optional fields
"""

OptionalFloat32: TypeAlias = Float32 | None
"""Optional 32-bit floating-point array type.

Extends Float32 to include None, representing memory-efficient arrays
that may be missing or not yet initialized.

Type Definition:
    `npt.NDArray[np.float32] | None`

See Also:
    Float32: Non-optional version
    OptionalFloat64: Higher precision alternative
"""

OptionalInt64: TypeAlias = Int64 | None
"""Optional 64-bit signed integer array type.

Extends Int64 to include None, representing large-range integer arrays
that may be missing or not yet initialized.

Type Definition:
    `npt.NDArray[np.int64] | None`

See Also:
    Int64: Non-optional version
    OptionalInt32: More memory-efficient alternative
"""

OptionalInt32: TypeAlias = Int32 | None
"""Optional 32-bit signed integer array type.

Extends Int32 to include None, representing integer arrays that may be
missing or not yet initialized.

Type Definition:
    `npt.NDArray[np.int32] | None`

Example:
    Optional connectivity data::

        bonds: Data[OptionalInt32] = Data[OptionalInt32](
            dtype=np.int32,
            meta=Metadata(description="Bond connectivity"),
            default=None
        )

See Also:
    Int32: Non-optional version
    OptionalInt64: Larger range alternative
"""

OptionalInt8: TypeAlias = Int8 | None

OptionalUInt64: TypeAlias = UInt64 | None

OptionalUInt32: TypeAlias = UInt32 | None

OptionalUInt8: TypeAlias = UInt8 | None
"""Optional 8-bit unsigned integer array type.

Extends Uint8 to include None, representing compact integer arrays
that may be missing or not yet initialized.

Type Definition:
    `npt.NDArray[np.uint8] | None`

Example:
    Optional atomic numbers::

        atom_z: Data[OptionalUint8] = Data[OptionalUint8](
            dtype=np.uint8,
            meta=Metadata(description="Atomic numbers"),
            default=None
        )

        # Check if atomic numbers are available
        z_values = molecule.atom_z.value
        if z_values is not None:
            unique_elements = np.unique(z_values)

See Also:
    Uint8: Non-optional version
    OptionalInt32: Signed alternative with larger range
"""

OptionalStr: TypeAlias = Str | None
"""Optional string array type.

Extends Str to include None, representing string arrays that may be
missing or not yet initialized.

Type Definition:
    `npt.NDArray[np.str_] | None`

Example:
    Optional element symbols::

        symbols: Data[OptionalStr] = Data[OptionalStr](
            dtype=np.str_,
            meta=Metadata(description="Element symbols"),
            default=None
        )

        # Conditional symbol processing
        symbol_data = molecule.symbols.value
        if symbol_data is not None:
            hydrogen_mask = symbol_data == 'H'

See Also:
    Str: Non-optional version
"""

DTYPE_TO_NUMPY: dict[type | TypeAlias, type] = {
    Float64: np.float64,
    Float32: np.float32,
    Int64: np.int64,
    Int32: np.int32,
    Int8: np.int8,
    UInt64: np.uint64,
    UInt32: np.uint32,
    UInt8: np.uint8,
    Str: np.str_,
    OptionalFloat64: np.float64,
    OptionalFloat32: np.float32,
    OptionalInt64: np.int64,
    OptionalInt32: np.int32,
    OptionalInt8: np.int8,
    OptionalUInt64: np.uint64,
    OptionalUInt32: np.uint32,
    OptionalUInt8: np.uint8,
    OptionalStr: np.str_,
}
