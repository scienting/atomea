# Expressions

In Atomea, at the heart of every selection lies an "expression."
An expression defines a specific criterion or a combination of criteria that atoms must satisfy to be included in a selection.
Think of expressions as the fundamental logical statements that evaluate to either `True` or `False` for each atom in your ensemble.

## What is an Expression?

Conceptually, an expression is a rule that operates on your `Ensemble` data to generate a boolean mask.
For example:

- "Is this atom an oxygen atom?"
- "Does this atom belong to molecule ID 5?"
- "Is this atom within 10 angstroms of a specific reference point?"

Each of these questions can be translated into an expression.
When an expression is "evaluated," it examines the relevant ensemble data (e.g., atom types, molecule IDs, coordinates) and produces a `True` or `False` value for every atom, indicating whether that atom meets the expression's condition.

## Ensemble Data

The power of expressions stems from their ability to interact directly with the various types of data stored within your `Ensemble`.
Different kinds of expressions are designed to leverage different "cadences" of data:

- **Ensemble-level data:** Some properties, like atom types or molecule IDs, typically remain constant across all microstates within a given run of an ensemble.
  Expressions based on these properties (e.g., `AtomTypeIs`, `MolIdIs`) can efficiently determine their truth value once per ensemble and apply it consistently.
- **Microstate-level data:** Other properties, most notably atomic coordinates, change from one microstate to another.
  Expressions that depend on these dynamic properties (e.g., `DistanceWithin`) must be evaluated independently for each microstate, yielding a unique mask for each snapshot in time.

This close coupling between expressions and the underlying ensemble data ensures that selections are accurate and reflect the true state of your system at any given moment or across the entire ensemble.

## Logical Operations

Just as you can combine simple logical statements in everyday language using "and," "or," and "not," Atomea allows you to combine individual expressions to form highly sophisticated selection queries.
This is achieved through logical operations:

- **AND (`&`)**: When two expressions are combined with an AND, an atom is selected only if it satisfies *both* expressions.
  For instance, "atoms that are oxygen AND are part of molecule ID 10."
- **OR (`|`)**: When two expressions are combined with an OR, an atom is selected if it satisfies *either* of the expressions (or both).
  For example, "atoms that are carbon OR are part of molecule ID 5."
- **NOT (`~`)**: The NOT operator inverts the result of an expression.
  If an expression would select an atom, applying NOT to it will cause that atom to be excluded, and vice-versa. For example, "atoms that are NOT oxygen."

These logical operations allow you to build up complex selection logic from simple, atomic expressions, providing immense flexibility in defining precisely the subset of data you need for your analysis.

