# Molecule IDs

Beyond individual atom properties, molecular simulations often require focusing on entire molecules or specific molecular groups.
Atomea facilitates this through selections based on molecule IDs, allowing you to easily isolate compounds, solvent molecules, or other predefined molecular entities within your ensemble.

The `MolIdIs` expression in Atomea enables you to select atoms that belong to molecules identified by a specific integer ID or a list of IDs.
Molecule IDs are typically assigned during system setup (e.g., in your topology file) and serve to group atoms that constitute a single molecule.
Like atom types, molecule IDs are generally static properties that remain constant across all microstates within a given ensemble.

When you utilize `MolIdIs`, the expression references the molecule ID information stored within your `Ensemble`.
It then generates a boolean mask where `True` indicates an atom belongs to one of the specified molecule IDs, and `False` indicates it does not.

