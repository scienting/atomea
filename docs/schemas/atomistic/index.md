# Atomistic Schemas

The atomistic schema in atomea is designed to provide a flexible, rigorous, and theory-driven interface for interacting with data originating from simulations and calculations on atomic systems.
In contrast to many domain-specific formats or narrowly scoped schemas (like those focused only on quantum chemistry, proteins, or molecular mechanics), atomea takes a pedantic and foundational approach.
Our schemas describe data by its theoretical or mathematical meaning, decoupled from the idiosyncrasies of file formats, domain assumptions, or tool-specific naming conventions.

This design provides interoperability, clarity, and precision, making it easier to build robust tools, standardize datasets, and bridge communities working across different subfields.

## Schema Design

Each ensemble is represented by an `EnsembleSchema`, which includes:

-   `identification`: metadata and fingerprinting;
-   `system`: atomic composition, charges, coordinates;
-   `topology`: atom types, bonding, and groupings;
-   `energy`: electronic, kinetic, and classical energies;
-   `qc`: quantum chemistry parameters and results;
-   `time`: time step and sampling intervals.

These fields are built from modular Pydantic models like `SystemSchema`, `EnergySchema`, `QCSchema`, and `TimeSchema`, all adhering to the same UUID/cadence principles.

Together, they form a rigorously structured yet highly extensible data model for representing atomic systems at scale.
