# Usage

Atomea's selection system provides a flexible and intuitive way to pinpoint specific atoms and molecules within your ensemble data. Whether you need to isolate a particular chemical species, focus on a binding site, or analyze interactions within a defined spatial region, selections are your go-to tool. This guide will walk you through the core concepts, available criteria, logical operators, and practical usage with detailed examples.1. Getting Started: The EnsembleSelectorAll selection operations in Atomea begin with the EnsembleSelector. This object acts as your query builder, allowing you to chain various criteria and logical operations to construct complex selection rules.To initialize an EnsembleSelector, you need to provide it with an Ensemble object and optionally a run_id if you wish to restrict your selection to a specific simulation run within that ensemble.

```python
from atomea.containers import Ensemble, Project
from atomea.selection.selector import EnsembleSelector

# Assuming you have an 'ensemble' object loaded from your project
# For demonstration, let's imagine 'my_ensemble' is already loaded.
# If you want to select from a specific run, provide its ID:
# selector = EnsembleSelector(my_ensemble, run_id="my_simulation_run_1")

# If you want to select across all runs (or if your ensemble only has one run):
selector = EnsembleSelector(my_ensemble)
```

Once initialized, the selector object becomes your canvas for building selection queries.

## Basic Selection Criteria

Atomea offers several fundamental methods to define your initial selection criteria. These methods return the EnsembleSelector object itself, allowing for seamless chaining.

### Selecting by Atom Type

This method allows you to select atoms based on their chemical element or a more specific atom type defined in your force field.
You provide a list of strings representing the desired atom types.

Example: Selecting all oxygen atoms (assuming "OW" is the atom type for water oxygen).

```python
# Select all oxygen atoms
oxygen_selector = EnsembleSelector(my_ensemble).atom_types(["OW"])

# You can also select multiple atom types
# Select all hydrogen atoms (assuming "HW" for water H, "H" for other H)
hydrogen_selector = EnsembleSelector(my_ensemble).atom_types(["HW", "H"])
```

### Selecting by Molecule ID

This method is useful for selecting entire molecules or groups of molecules based on their assigned integer IDs.
Example: Selecting all atoms belonging to molecule ID 0 (often the first molecule in a system).

```python
# Select all atoms in molecule with ID 0
mol0_selector = EnsembleSelector(my_ensemble).molecule_ids([0])

# Select atoms in molecules with IDs 0 and 2
mols_0_and_2_selector = EnsembleSelector(my_ensemble).molecule_ids([0, 2])
```

### Selecting by Distance

This powerful method allows you to select atoms that are within a specified radial distance from a reference point or a group of reference atoms.
The distance is measured in the same units as your ensemble's coordinates (e.g., Angstroms).
The from_atoms argument can be specified in a few ways:

- A single atom index (integer): Selects atoms within dist of a specific atom.
- A slice of atom indices (slice object): Selects atoms within dist of any atom within that index range.
- Another EnsembleSelector instance: This is incredibly flexible, allowing you to define your reference group using any other selection criteria.

Atoms will be selected if they are within dist of any atom in the reference group.

!!! important
    Important Note on Cadence: Unlike atom types and molecule IDs, atomic coordinates change with each microstate. Therefore, distance_within selections are evaluated per microstate, yielding a unique mask for each snapshot in your ensemble.

Examples:

```python
# Select atoms within 5.0 Angstroms of atom with index 15
near_atom_15_selector = EnsembleSelector(my_ensemble).distance_within(from_atoms=15, dist=5.0)

# Select atoms within 3.0 Angstroms of any atom from index 0 to 9 (inclusive)
near_first_10_atoms_selector = EnsembleSelector(my_ensemble).distance_within(from_atoms=slice(0, 10), dist=3.0)

# Select atoms within 4.0 Angstroms of any atom that is a "protein" atom
# (Assuming you have a way to define "protein" atoms, e.g., by molecule ID or residue type)
# First, define the reference selection:
protein_atoms_ref = EnsembleSelector(my_ensemble).molecule_ids([1, 2, 3]) # Example: Mol IDs 1, 2, 3 are protein
near_protein_selector = EnsembleSelector(my_ensemble).distance_within(from_atoms=protein_atoms_ref, dist=4.0)
```

## Combining Selections

The true power of Atomea's selection system comes from its ability to combine basic criteria using logical operators: AND, OR, and NOT.
These operators allow you to build highly specific and complex queries.

### Implicit AND (Chaining)

When you chain multiple selection methods directly, they are implicitly combined with a logical AND.
This means an atom must satisfy all chained conditions to be selected.
Example: Select carbon atoms that are also part of molecule ID 2.

```python
# Select atoms that are 'C' type AND belong to molecule ID 2
carbon_in_mol2_selector = EnsembleSelector(my_ensemble).atom_types(["C"]).molecule_ids([2])
# This is equivalent to: (atom_type is 'C') AND (molecule_id is 2)
```

### Explicit AND()

While chaining implies an AND, you can explicitly use the .AND() method for clarity, especially in more complex expressions or when you want to visually separate logical blocks.
It does not change the behavior of implicit AND.

Example: Explicitly showing the AND operation.

```python
# Select atoms that are 'C' type AND belong to molecule ID 2 (explicit AND)
carbon_in_mol2_explicit_selector = (
    EnsembleSelector(my_ensemble)
    .atom_types(["C"])
    .AND() # Explicitly state the AND
    .molecule_ids([2])
)
```

### Logical OR

The .OR() method combines the preceding selection criteria with the next one using a logical OR.
An atom is selected if it satisfies either the accumulated criteria or the next criterion.

!!! important
  You cannot start a selection query with .OR(). You must have at least one selection criterion defined before using OR().


Example: Select all oxygen atoms OR all atoms belonging to molecule ID 2.

```python
# Select atoms that are 'OW' type OR belong to molecule ID 2
oxygen_or_mol2_selector = (
    EnsembleSelector(my_ensemble)
    .atom_types(["OW"])
    .OR()
    .molecule_ids([2])
)
# This is equivalent to: (atom_type is 'OW') OR (molecule_id is 2)
```

### Logical NOT

The .NOT() method applies a logical NOT operation to the next selection criterion that follows it.
This means it selects atoms that do not satisfy that criterion.

```python
Example: Select all atoms that are NOT carbon.# Select all atoms that are NOT 'C' type
not_carbon_selector = EnsembleSelector(my_ensemble).NOT().atom_types(["C"])

# Select atoms in molecule 0 AND (NOT oxygen atoms)
mol0_and_not_oxygen_selector = (
    EnsembleSelector(my_ensemble)
    .molecule_ids([0])
    .AND()
    .NOT()
    .atom_types(["OW"])
)
# This is equivalent to: (molecule_id is 0) AND (NOT (atom_type is 'OW'))
```

## Retrieving the Selection

After you've built your selection query using the EnsembleSelector and its methods, the final step is to execute it and retrieve the results.
This is done by calling the .get_mask() method.
get_mask() performs the lazy evaluation of your selection expression and yields an iterator of boolean NumPy arrays.
Each array represents a "mask" for a single microstate within your ensemble:

- True at an atom's position means that atom is included in the selection for that microstate.
- False means it is excluded.

You can control which microstates you get masks for using the micro_id argument:

- None (default): Yields masks for all microstates in the ensemble.
- An integer: Yields the mask for a single, specific microstate (e.g., micro_id=0 for the first microstate).
A slice object: Yields masks for a range of microstates (e.g., micro_id=slice(10, 20) for microstates 10 through 19).

Example: Getting masks for a selection.

```python
import numpy as np

# Let's use a selection that finds all oxygen atoms
oxygen_selector = EnsembleSelector(my_ensemble).atom_types(["OW"])

# Get masks for all microstates
all_oxygen_masks = oxygen_selector.get_mask(micro_id=slice(None))

print("Masks for all microstates:")
for i, mask in enumerate(all_oxygen_masks):
    print(f"  Microstate {i}: {mask}")
    # You can then use this mask for further analysis or visualization
    # For example, to get the indices of selected atoms:
    # selected_indices = np.where(mask)[0]
    # print(f"    Selected atom indices: {selected_indices}")

# Get mask for a specific microstate (e.g., the first one)
first_microstate_oxygen_mask = oxygen_selector.get_mask(micro_id=0)

print("\nMask for the first microstate:")
for mask in first_microstate_oxygen_mask: # This will yield only one mask
    print(f"  Microstate 0: {mask}")
```

## Putting It All Together

Here are some more elaborate examples demonstrating how to combine different criteria and logical operators to achieve precise selections.

### Water molecules near a specific carbon atom.

Select all atoms belonging to water molecules (Mol IDs 0 and 1, assuming) that are also within 6.0 Angstroms of atom 6 (a carbon atom).

```python
# Define the reference carbon atom (atom with index 6)
carbon_ref = 6

# Build the complex selector
water_near_carbon_selector = (
    EnsembleSelector(my_ensemble)
    .molecule_ids([0, 1]) # Selects atoms in water molecules
    .AND()
    .distance_within(from_atoms=carbon_ref, dist=6.0) # And are near the carbon
)

print("\nWater atoms near carbon (complex selection):")
for i, mask in enumerate(water_near_carbon_selector.get_mask()):
    print(f"  Microstate {i}: {mask}")
```

### Atoms that are either Oxygen OR (Hydrogen AND in Molecule 2)

This demonstrates the use of parentheses for grouping logical operations, which is achieved implicitly through chaining and explicit OR/AND calls.

```python
# Select Oxygen atoms
oxygen_atoms = EnsembleSelector(my_ensemble).atom_types(["OW"])

# Select Hydrogen atoms that are in Molecule 2
hydrogen_in_mol2 = EnsembleSelector(my_ensemble).atom_types(["H", "HW"]).molecule_ids([2])

# Combine them with OR: (Oxygen) OR (Hydrogen AND in Molecule 2)
combined_selector = oxygen_atoms.OR()._add_expression(hydrogen_in_mol2._current_expression) # Accessing internal expression for complex OR

# A more direct way to achieve this with the fluent API would be:
# (Note: The current fluent API design for OR chaining needs careful consideration
# when combining already built sub-selectors. For simple cases, direct chaining works.)
# For this specific case, if `hydrogen_in_mol2` was built first, you'd do:
# combined_selector = (
#     EnsembleSelector(my_ensemble)
#     .atom_types(["OW"])
#     .OR()
#     .atom_types(["H", "HW"])
#     .AND()
#     .molecule_ids([2])
# )
# This would evaluate to: (OW) OR (H/HW AND Mol2)

# Let's use the more explicit way for clarity of the concept:
# Imagine `sub_selector_A` is `atom_types(["OW"])`
# Imagine `sub_selector_B` is `atom_types(["H", "HW"]).molecule_ids([2])`
# Then you'd want `sub_selector_A.OR()._add_expression(sub_selector_B._current_expression)`
# For simplicity, let's assume a direct fluent chain for this example as the API is designed:
complex_selector_fluent = (
    EnsembleSelector(my_ensemble)
    .atom_types(["OW"]) # Start with Oxygen
    .OR()                # OR with the next group
    .atom_types(["H"])   # This starts the second group (Hydrogen)
    .molecule_ids([2])   # Which is ANDed with Mol ID 2
)
# This reads as: (atom_type is 'OW') OR ((atom_type is 'H') AND (molecule_id is 2))

print("\nComplex selector (Oxygen OR (Hydrogen AND Mol 2)):")
for i, mask in enumerate(complex_selector_fluent.get_mask()):
    print(f"  Microstate {i}: {mask}")
```

