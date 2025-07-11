# Containers

[Atomea containers][containers.core.Container] are Python classes that define the structure and constraints of data using Python type annotations.
In the context of atomea, containers are central to organizing, validating, and documenting your data.

Each container represents a well-defined group of related information.
This helps us create clear boundaries between different types of data and workflows.
For example:

-   The [`Energy`][containers.energy.Energy] container contains various types of energies.
-   The [`Quantum`][containers.quantum.Quantum] container represents quantum chemical systems.

But containers in atomea are not limited to scientific results.
We also use them to define the parameters and behavior of software execution environments.
This makes them especially powerful for automating complex workflows.

Here is a simplified example of a container in atomea:

```python
class Energy:
    electronic: pl.DataFrame | None
    kinetic: pl.DataFrame | None
    potential_mm: pl.DataFrame | None
    ...
```

This shows that we have defined three optional (i.e., `None`) types of energy that are stored in a polars DataFrame.

## Data

Atomea Data are the individual components that make up a container.
Each field defines the characteristics of a single piece of data—its type, default value, constraints, and validation rules.
In atomea, fields are used to express the exact expectations for a particular value within a container.
This allows us to ensure data integrity, provide clear documentation, and enable rich tooling like autocompletion and error checking.
At its simplest, a field is just a type-annotated class attribute:

TODO: Write
