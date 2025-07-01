# Containers

Atomea containers are Python classes that define the structure and constraints of data using Python type annotations.
In the context of atomea, containers are central to organizing, validating, and documenting scientific and computational data.

Each container represents a well-defined group of related information.
This helps us create clear boundaries between different types of data and workflows.
For example:

-   The `Energy` container contains various types of energies (e.g., total energy, kinetic energy, zero-point energy).
-   The `Quantum` container represents quantum chemical systems, including their molecular geometry, basis sets, and computational methods.

But containers in atomea aren't limited to scientific results.
We also use them to define the parameters and behavior of software execution environments.
This makes them especially powerful for automating complex workflows.
For instance: to run molecular dynamics simulations with Amber, you need to specify many input parameters.
Remembering all the defaults and syntax can be tedious.
With atomea, these inputs are encapsulated in a container that serves both as documentation and as a validation tool, ensuring your simulation configuration is complete and correct.

Here's a simplified example of a container in atomea:

```python
class EnergySchema:
    total_energy: float
    kinetic_energy: float | None = None
    zero_point_energy: float | None = None
```

This guarantees that `total_energy` must always be present and must be a float, while the other two fields are optional.

## Data

Atomea Data are the individual components that make up a container.
Each field defines the characteristics of a single piece of data—its type, default value, constraints, and validation rules.
In atomea, fields are used to express the exact expectations for a particular value within a container.
This allows us to ensure data integrity, provide clear documentation, and enable rich tooling like autocompletion and error checking.
At its simplest, a field is just a type-annotated class attribute:

TODO: Write
