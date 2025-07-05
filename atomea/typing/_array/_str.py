from typing import TypeAlias

import numpy as np
import numpy.typing as npt

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
        symbol_data = molecule.symbols.values()
        if symbol_data is not None:
            hydrogen_mask = symbol_data == 'H'

See Also:
    Str: Non-optional version
"""
