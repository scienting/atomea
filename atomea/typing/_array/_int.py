from typing import TypeAlias

import numpy as np
import numpy.typing as npt

from atomea.data import Data

Int64NP: TypeAlias = npt.NDArray[np.int64]

Int64: TypeAlias = Data[Int64NP]
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


Int32NP: TypeAlias = npt.NDArray[np.int32]


Int32: TypeAlias = Data[Int32NP]
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

Int16NP: TypeAlias = npt.NDArray[np.int16]

Int16: TypeAlias = Data[Int16NP]


Int8NP: TypeAlias = npt.NDArray[np.int8]

Int8: TypeAlias = Data[Int8NP]

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

OptionalInt16: TypeAlias = Int16 | None

OptionalInt8: TypeAlias = Int8 | None
