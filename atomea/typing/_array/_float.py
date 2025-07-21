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
        coord_data: Float64 = ensemble.coordinates.read()

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

        energy_data = ensemble.energy.read()
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
