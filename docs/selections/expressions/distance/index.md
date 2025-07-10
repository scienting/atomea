# Distance

Spatial proximity is a critical aspect of analyzing molecular systems, allowing you to investigate interactions, define active sites, or identify solvent shells around solutes. Atomea's "distance within" selections provide a robust way to identify atoms located within a specified radius of a reference point or a set of other atoms.

The `DistanceWithin` expression in Atomea allows you to select atoms that are located within a certain `dist` (distance) from a specified `from_atoms` reference. This reference can be:

- A single atom identified by its index.
- A slice of atoms, specifying a range of atom indices.
- The result of another `SelectionExpression`, meaning you can define your reference group using any other Atomea selection criterion (e.g., all "protein" atoms, or specific ligand atoms).

Unlike atom types or molecule IDs, atomic coordinates are dynamic, changing from one microstate (snapshot) to the next in your ensemble. Consequently, a `DistanceWithin` expression is evaluated independently for each microstate. This ensures that the selection accurately reflects the spatial arrangement of atoms at every point in time, producing a unique boolean mask for each microstate. The mask indicates `True` for atoms that fall within the specified distance from *any* of the reference atoms in that particular microstate, and `False` otherwise.

