# Microstate

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
class Energy:
    electronic: ...
```

An ensemble can have zero or many microstates.
