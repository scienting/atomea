from typing import TypeAlias

import polars as pl

DataFrame: TypeAlias = pl.DataFrame
"""Deprecated: Use OptionalDataFrame instead.

This alias is maintained for backward compatibility but OptionalDataFrame
should be used for new code to be explicit about the optional nature.

Type Definition:
    `pl.DataFrame | None`

Warning:
    This alias may be removed in future versions. Use OptionalDataFrame
    for new code.

See Also:
    OptionalDataFrame: Preferred alternative
"""

OptionalDataFrame: TypeAlias = DataFrame | None
"""Optional Polars DataFrame type.

Represents tabular data that may be missing or not yet initialized.
DataFrames are used for structured data that doesn't fit naturally
into multidimensional arrays, such as analysis results, metadata
tables, or relational data.

Type Definition:
    `pl.DataFrame | None`

Typical Use Cases:
    - Analysis results with multiple columns
    - Metadata tables with mixed data types
    - Time series data with irregular intervals
    - Aggregated statistics and summaries
    - Relational data requiring joins and grouping

Example:
    Optional analysis results::

        analysis: Data[OptionalDataFrame] = Data[OptionalDataFrame](
            dtype=None,  # DataFrames don't have numpy dtypes
            meta=Metadata(
                description="Statistical analysis results",
                store=StoreKind.TABLE
            ),
            default=None
        )

    Creating DataFrame data::

        import polars as pl

        results = pl.DataFrame({
            'atom_id': [1, 2, 3],
            'element': ['H', 'H', 'O'],
            'charge': [-0.4, -0.4, -0.8],
            'energy': [1.2, 1.3, -5.4]
        })
        ensemble.analysis = results

    Conditional processing::

        data = ensemble.analysis.read()
        if data is not None:
            hydrogen_data = data.filter(pl.col('element') == 'H')
            mean_energy = data['energy'].mean()

DataFrame Advantages:
    - Mixed data types in columns
    - Efficient columnar storage
    - Rich query and transformation API
    - Memory efficient for sparse data
    - SQL-like operations

See Also:
    DataFrame: Non-optional alias (deprecated, use this instead)
    polars.DataFrame: Polars DataFrame documentation
    atomea.stores.TableStore: DataFrame storage backend
"""
