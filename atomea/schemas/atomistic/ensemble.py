from typing import Any, Generator

from pydantic import BaseModel, Field

from ..data import Data
from ..id import IdentificationSchema
from ..io import YamlIO
from .qc import QCSchema
from .system import SystemSchema
from .topology import TopologySchema


class EnsembleSchema(BaseModel, YamlIO, Data):
    """
    The `EnsembleSchema` class is a Pydantic model designed to represent a
    collection of molecular structures, referred to as frames. This class is used to
    manage and validate an ensemble of molecular data, facilitating the handling of
    multiple molecular configurations, such as those produced during atomistic
    calculations.
    """

    identification: IdentificationSchema = Field(default_factory=IdentificationSchema)
    """
    This attribute stores identification information for the molecule, such as its
    name, unique identifiers, and other metadata. It uses the IdentificationSchema
    to validate and structure the data.
    """

    qc: QCSchema = Field(default_factory=QCSchema)
    """
    This attribute holds quantum chemistry (QC) properties of the molecule, such as
    energy, wavefunction information, and other relevant QC data. It is structured and
    validated using the QCSchema.
    """

    system: SystemSchema = Field(default_factory=SystemSchema)
    """
    This attribute contains system-related properties of the molecule, such as
    atomic coordinates, velocities, and other system-specific data. The SystemSchema
    is used for validation and structuring of this data.
    """

    topology: TopologySchema = Field(default_factory=TopologySchema)
    """
    This attribute includes topological information about the molecule, such as the
    bonding structure, atom types, and other topology-related data. The TopologySchema
    validates and structures this data.
    """

    @classmethod
    def generate_fields(
        cls, model: Any, parent_key: str = ""
    ) -> Generator[tuple[str, Any], None, None]:
        """
        Recursively generates all FieldInfo objects from a nested Pydantic BaseModel.

        Args:
            model: The Pydantic model instance to process.

        Yields:
            A generator yielding ModelField objects.
        """
        for field_name, field in model.__fields__.items():
            key = f"{parent_key}.{field_name}" if parent_key else field_name
            yield key, field
            field_value = getattr(model, field_name)
            if isinstance(field_value, BaseModel):
                yield from cls.generate_fields(field_value, key)

    @classmethod
    def get_schema_map(cls, model: BaseModel) -> dict[str, dict[str, str]]:
        """
        Generates a mapping of full field keys to UUIDs from a nested Pydantic
        BaseModel.

        Args:
            model: The Pydantic model instance to process.

        Returns:
            Keys are the UUID of a single field and values are

                -   `field_key`: key to the model field from
                    [`generate_fields`]
                    [schemas.atomistic.ensemble.EnsembleSchema.generate_fields].
                -   `cadence`: The frequency of when this data would change; per molecule
                    or ensemble basis.
        """
        uuid_mapping = {}
        for key, field in cls.generate_fields(model):
            if not hasattr(field, "metadata") or len(field.metadata) == 0:
                continue
            uuid = field.metadata[0].get("uuid", None)
            cadence = field.metadata[0].get("cadence", None)
            if uuid:
                uuid_mapping[uuid] = {"field_key": key, "cadence": cadence}
        return uuid_mapping
