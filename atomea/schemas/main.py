from typing import Any

from pydantic import BaseModel, Field

from .id import IdentificationSchema
from .qc import QCSchema
from .system import SystemSchema
from .topology import TopologySchema


class MoleculeSchema(BaseModel):
    """The MoleculeSchema class is a Pydantic model designed to represent a molecule in
    a structured and validated format. This class integrates various aspects of a
    molecule, including its identification, quantum chemistry properties,
    system properties, and topology. By leveraging Pydantic's data validation
    capabilities, MoleculeSchema ensures that the data adheres to the expected
    structure and types, facilitating reliable and consistent data handling in
    computational chemistry and biology workflows.

    The MoleculeSchema class is essential for managing and validating complex molecular
    data in computational chemistry and biology projects. By providing a structured and
    validated representation of a molecule, it facilitates various tasks such as:

    -   Data integration from multiple sources.
    -   Consistent data handling and processing.
    -   Simplified data exchange between different components of a computational pipeline.
    -   Enhanced data validation to ensure the integrity and correctness of molecular data.

    The update method enhances the flexibility of the MoleculeSchema class by allowing
    dynamic updates to its attributes, making it adaptable to changing data requirements
    and facilitating seamless data manipulation.
    """

    identification: IdentificationSchema = Field(default_factory=IdentificationSchema)
    """This attribute stores identification information for the molecule, such as its
    name, unique identifiers, and other metadata. It uses the IdentificationSchema
    to validate and structure the data.
    """

    qc: QCSchema = Field(default_factory=QCSchema)
    """This attribute holds quantum chemistry (QC) properties of the molecule, such as
    energy, wavefunction information, and other relevant QC data. It is structured and
    validated using the QCSchema."""

    system: SystemSchema = Field(default_factory=SystemSchema)
    """This attribute contains system-related properties of the molecule, such as
    atomic coordinates, velocities, and other system-specific data. The SystemSchema
    is used for validation and structuring of this data.
    """

    topology: TopologySchema = Field(default_factory=TopologySchema)
    """This attribute includes topological information about the molecule, such as the
    bonding structure, atom types, and other topology-related data. The TopologySchema
    validates and structures this data."""

    def update(self, data: dict[str, Any]) -> None:
        """Update the fields of the MoleculeSchema instance with the provided data.

        This method updates the attributes of the MoleculeSchema instance based on
        the keys and values in the provided dictionary. The keys in the dictionary
        can represent nested fields using dot notation.

        Args:
            data (dict[str, Any]): A dictionary containing the keys and values to
            update the MoleculeSchema instance. The keys can use dot notation to
            specify nested attributes.

        Example:
            >>> molecule = MoleculeSchema()
            >>> update_data = {
            >>>     "identification.name": "Water",
            >>>     "qc.energy": -76.4,
            >>>     "system.coordinates": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            >>>     "topology.bonds": [(0, 1), (0, 2)]
            >>> }
            >>> molecule.update(update_data)

        Raises:
            AttributeError: If a specified attribute does not exist in the schema.

        Notes:
            -   The method supports updating nested attributes by splitting keys on the
                dot ('.') character.
            -   If a key does not use dot notation, it will update the top-level
                attribute directly.
        """
        for key, value in data.items():
            keys = key.split(".")
            if len(keys) > 1:
                sub_model = self
                for sub_key in keys[:-1]:
                    sub_model = getattr(sub_model, sub_key)
                setattr(sub_model, keys[-1], value)
            else:
                setattr(self, key, value)


class EnsembleSchema(BaseModel):
    """The EnsembleSchema class is a Pydantic model designed to represent a
    collection of molecular structures, referred to as frames. This class is used to
    manage and validate an ensemble of molecular data, facilitating the handling of
    multiple molecular configurations, such as those produced during molecular dynamics
    simulations or geometry optimizations.
    """

    frames: list[MoleculeSchema] = Field(default=[])
    """The frames attribute is a list that stores instances of MoleculeSchema. Each
    instance represents a single molecular structure or configuration within the
    ensemble. This attribute allows the EnsembleSchema to manage multiple molecular
    structures collectively, making it easier to handle data from simulations that
    produce multiple frames, such as molecular dynamics trajectories or conformational
    scans."""
