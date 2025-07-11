# Arrays

As discussed [before](/data/dimensionality), data can have more than one dimensions.
We use an [`ArrayStore`][stores.arrays.core.ArrayStore] when dealing with multi-dimensional data (e.g., [`ZarrArrayStore`][stores.arrays._zarr.ZarrArrayStore] and [`NumpyArrayStore`][stores.arrays._numpy.NumpyArrayStore]) to ensure high performance with numerical data.

## Backends

Remember, atomea is designed so that [users don't have to worry about how data is stored](/stores/); thus, interacting with atomea data is the exact same for every backend.
Choosing a backend is made based on the specific use cases and all data can be seamlessly converted between formats.

### Zarr

TODO: Introduce the who and where of zarr.
Discuss when to or not to choose this array backend.

### NumPy

TODO: Introduce the who and where of NPY and NPZ.
Discuss when to or not to choose this array backend.
