# Data by Meaning, Not by Name

In scientific practice, the term "energy" may refer to many things such as electronic energy, Gibbs free energy, or molecular mechanics potential energy.
Many tools and datasets use the same word to mean different things—or different words to mean the same thing.
This inconsistency becomes a major barrier to integration.

Atomea resolves this by explicitly defining each piece of data using its theoretical definition.
Every field is:

-   Described precisely (e.g., electronic energy as the expectation value of the electronic Hamiltonian);
-   Annotated with a cadence, specifying at what level the data applies (e.g., per `ensemble` or per `microstate`).

Every field in the schema is also assigned a universally unique identifier (UUID) that looks like this `b76e16e6-3d22-4e0a-80c9-9f0901dc8c64`.
This serves multiple critical purposes:

-   Programs can communicate unambiguously without relying on field names.
-   Users can alias field names freely (e.g., rename `"electronic"` to `"E_el"`), and atomea will still know exactly what type of data it refers to.
-   You can query and transform datasets by theoretical quantity, not just by label.

UUIDs also power atomea’s ability to interoperate with external datasets: as long as a dataset includes the UUID, it can be validated, merged, or interpreted correctly regardless of its original format.

