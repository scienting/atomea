from pydantic import BaseModel, Field

from ..id import IdentificationSchema
from ..io import YamlIO
from .qc import QCSchema
from .system import SystemSchema
from .topology import TopologySchema


class MoleculeSchema(BaseModel, YamlIO):
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
