# Atom type

Identifying atoms based on their chemical type (e.g., carbon, oxygen, hydrogen) is a fundamental and frequently performed operation.
Atomea provides a straightforward way to achieve this through dedicated selection expressions.

The `AtomTypeIs` expression in Atomea allows you to select atoms whose chemical type matches one or more specified types.
Atom types are typically defined in your system's topology and represent the fundamental elemental or force-field-specific character of each atom.
This property is generally static across all microstates within an ensemble, meaning an atom's type does not change over the course of a simulation.

When you use `AtomTypeIs`, the expression examines the atom type information stored within your `Ensemble` and generates a boolean mask. For each atom, this mask indicates `True` if its type is found within the list you provided, and `False` otherwise.

