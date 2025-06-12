# Schema Design

Each ensemble is represented by an `Ensemble`, which includes:

-   `identification`: metadata and fingerprinting;
-   `system`: atomic composition, charges, coordinates;
-   `topology`: atom types, bonding, and groupings;
-   `energy`: electronic, kinetic, and classical energies;
-   `qc`: quantum chemistry parameters and results;
-   `time`: time step and sampling intervals.

These fields are built from modular Pydantic models like `Microstate`, `EnergySchema`, `QCSchema`, and `TimeSchema`, all adhering to the same UUID/cadence principles.
Together, they form a rigorously structured yet highly extensible data model for representing atomic systems at scale.

::: containers.ensemble.Ensemble
    handler: python
    options:
      inherited_members: true
      show_root_heading: true
      show_root_full_path: false
      show_root_members_full_path: false
      show_labels: false
