# Ensemble

An ensemble defines a molecular system and all the conditions that remain fixed for its entirety:

-   Force field parameters;
-   Simulation temperature;
-   Molecular topology;
-   Environmental assumptions (e.g., solvent, pH);
-   Scoring functions for docking, etc.

These are defined by fields that use `cadence: ensemble`, and are meant to remain constant within that ensemble.

A new ensemble must be created if any of these assumptions change.
For example:

-   A molecular dynamics simulation run at pH 7 is one ensemble.
-   The same system at pH 5 is a different ensemble.
