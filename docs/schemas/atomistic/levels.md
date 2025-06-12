# Ensembles and Microstates

Borrowing from statistical mechanics, we organize atomistic data across two key levels.

## Ensemble

An ensemble defines a molecular system and all the conditions that remain fixed for its entirety:

-   Force field parameters;
-   Simulation temperature;
-   Molecular topology;
-   Environmental assumptions (e.g., solvent, pH);
-   Scoring functions for docking, etc.

These are defined by fields that use `cadence: ensemble`, and are meant to remain constant within that ensemble.

For example:

```python
from typing import Annotated

class Microstate(BaseModel):
    atom_z: Annotated[np.ndarray, {"cadence": "ensemble", "uuid": "d051abd9-c815-40b1-ab2d-e7a50a2d3259"}]
    ...
```

A new ensemble must be created if any of these assumptions change.
For example:

-   A molecular dynamics simulation run at pH 7 is one ensemble.
-   The same system at pH 5 is a different ensemble.

## Microstate

A microstate is a single atomistic configuration—a snapshot in space and time:

-   One frame in a trajectory;
-   One pose from a docking simulation;
-   One QM calculation of a conformation.

Fields at this level use `cadence: microstate` and describe things that can vary across microstates:

-   Coordinates;
-   Velocities;
-   Energies;
-   Electronic structure.

For instance:

```python
from typing import Annotated

class EnergySchema(BaseModel):
    electronic: Annotated[np.ndarray | None, {"cadence": "microstate", "uuid": "9e4bdf45-0150-4605-9528-e23aed0be9f2"}]
```

An ensemble can have zero or many microstates.

## Examples of Use Cases

-   Molecular Dynamics: A 1000-frame MD trajectory produces one ensemble with 1000 microstates.
-   Docking Study: Docking 10 ligands and keeping the best pose results in 10 ensembles with 1 microstate each.
-   Quantum Chemistry Benchmarking: Comparing different functionals on the same molecule means multiple ensembles (one per method), each with the same coordinates but different computed properties.


