# Dimensionality

The way we store and analyze data is significantly impacted by its dimensionality.
Think of dimensionality as the number of "directions" or "axes" needed to describe a piece of data.
Data can generally be categorized as either scalars or arrays.

### Scalars

When data is represented by a single value, we call it a scalar.
This means it has only one dimension.
Examples include things like a material's energy, a specific thermodynamic variable, or a calculation parameter.
Each of these can be fully described by just one number.

### Arrays

Data that consists of more than one value needs to be stored using arrays, such as NumPy arrays.
These arrays are crucial because they can hold multiple values and are structured with the appropriate number of dimensions to reflect that.
Even if you only have a few related values, they're still best stored in an array.
