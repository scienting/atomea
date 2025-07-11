# Tables

Much of our data are represented as scalars such as energies and scoring functions.
While these data can also benefit from being stored in [arrays](/stores/arrays/), we often want to be able to easily query the data across many [ensembles](/data/cadence#ensemble) or [runs](/data/cadence/#run).
Tables are inherently optimized for these types of queries with the added benefit of providing multiple columns as well.
