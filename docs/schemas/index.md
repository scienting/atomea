# Schemas

Schemas are a core concept for atomea and any other package that deals with data.
Our schemas are implemented using [Pydantic][pydantic] and serve as a way to document, validate, and interact with data.
In particular, our focus is on

-   [Atomistic schemas](./atomistic) that describe the data originating from various calculations of atomistic systems.
-   [Workflow schemas](./workflow) that assist with setting up computational chemistry and biology calculations.

Before we dive in, it is important to understand the core components of [Pydantic][pydantic] schemas: [models][pydantic-models] and [fields][pydantic-fields].

## Models

Pydantic [models][pydantic-models] are Python classes that define the structure and constraints of data using Python type annotations.
In the context of atomea, models are central to organizing, validating, and documenting scientific and computational data.

Each model represents a well-defined group of related information.
This helps us create clear boundaries between different types of data and workflows.
For example:

-   The `EnergySchema` model contains various types of energies (e.g., total energy, kinetic energy, zero-point energy).
-   The `QCSchema` model represents quantum chemical systems, including their molecular geometry, basis sets, and computational methods.

But models in atomea aren't limited to scientific results.
We also use them to define the parameters and behavior of software execution environments.
This makes them especially powerful for automating complex workflows.
For instance:

-   To run molecular dynamics simulations with Amber, you need to specify many input parameters.
    Remembering all the defaults and syntax can be tedious.
    With atomea, these inputs are encapsulated in a Pydantic model that serves both as documentation and as a validation tool, ensuring your simulation configuration is complete and correct.
-   Similarly, to run jobs on high-performance computing (HPC) clusters using Slurm, we provide a `SlurmSchema` model.
    This model captures the structure of a Slurm job submission script.
    You can dynamically define and update this schema in Python, and atomea can then render it into a ready-to-run Slurm script.

Here's a simplified example of a model in atomea:

```python
from pydantic import BaseModel

class EnergySchema(BaseModel):
    total_energy: float
    kinetic_energy: float | None = None
    zero_point_energy: float | None = None
```

This guarantees that `total_energy` must always be present and must be a float, while the other two fields are optional. Pydantic automatically validates the types, fills in defaults, and provides helpful error messages when input data is malformed.

## Fields

Pydantic [fields][pydantic-fields] are the individual components that make up a model.
Each field defines the characteristics of a single piece of data—its type, default value, constraints, and validation rules.
In atomea, fields are used to express the exact expectations for a particular value within a model.
This allows us to ensure data integrity, provide clear documentation, and enable rich tooling like autocompletion and error checking.
At its simplest, a field is just a type-annotated class attribute:

```python
from pydantic import BaseModel

class MoleculeSchema(BaseModel):
    name: str
    charge: int = 0
```

In this example:

-   `name` is a required field of type `str`.
-   `charge` is an optional field of type `int` with a default value of `0`.

But fields can do much more than that. Using the `Field()` function, you can add:

-   Validation constraints (e.g., min/max values, regex patterns)
-   Metadata (e.g., descriptions, units)
-   Examples or aliases for documentation and flexibility

```python
from pydantic import BaseModel, Field

class SimulationSettings(BaseModel):
    temperature: float = Field(..., gt=0, description="Temperature in Kelvin")
    timestep: float = Field(1.0, ge=0.1, le=5.0, description="Time step in femtoseconds")
```

Here:

-  `temperature` must be a positive float and is required (`...` indicates no default).
-  `timestep` must be between `0.1` and `5.0` fs, with a default of `1.0`.

Fields can also incorporate custom validators when standard constraints aren’t sufficient.
This is helpful for domain-specific rules—for example, ensuring that a string field matches a specific format used in quantum chemistry input files.

<!-- REFERENCES -->

[pydantic]: https://docs.pydantic.dev/latest/
[pydantic-models]: https://docs.pydantic.dev/latest/concepts/models/
[pydantic-fields]: https://docs.pydantic.dev/latest/concepts/fields/
