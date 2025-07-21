# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import TypeAlias

import numpy as np
import numpy.typing as npt

UInt64: TypeAlias = npt.NDArray[np.uint64]

UInt32: TypeAlias = npt.NDArray[np.uint32]

UInt16: TypeAlias = npt.NDArray[np.uint16]

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

OptionalUInt64: TypeAlias = UInt64 | None

OptionalUInt32: TypeAlias = UInt32 | None

OptionalUInt16: TypeAlias = UInt16 | None

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
        z_values = molecule.atom_z.read()
        if z_values is not None:
            unique_elements = np.unique(z_values)

See Also:
    Uint8: Non-optional version
    OptionalInt32: Signed alternative with larger range
"""
