from typing import Any, Collection, Generator

from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor

from loguru import logger

from ..schemas.atomistic import EnsembleSchema, MoleculeSchema


class Digester(ABC):
    """The Digester class is an abstract base class designed to assist with digesting
    data from various computational chemistry and biology resources. This class
    provides a framework for processing, parsing, and validating data extracted from
    simulations, geometry optimizations, and other computational methods.
    """

    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def checks(cls):
        """
        Perform basic checks to raise warnings or errors before digesting.

        This method should be overridden in the child class to include specific
        checks required for the particular type of data being processed.
        """
        pass

    @classmethod
    @abstractmethod
    def get_uuid_map(cls):
        """
        Update the function UUID map by inspecting the class methods
        decorated with @SchemaUUID.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def prepare_inputs_digester(
        cls, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        """Prepare and return the inputs necessary to start the digestion process.

        This abstract method must be implemented in any child class. It should
        return a dictionary of inputs that will be used by
        [`get_inputs_frame`][digesters.digester.Digester.get_inputs_frame].

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            A dictionary of inputs for the frame digesting process.
        """
        raise NotImplementedError

    @classmethod
    def digest(
        cls,
        digester_args: tuple[Any, ...] = tuple(),
        digester_kwargs: dict[str, Any] = dict(),
        parallelize: bool = False,
        max_workers: int | None = None,
    ) -> EnsembleSchema:
        """
        Given the inputs in `digester_args` and `digester_kwargs`, digest all possible
        atomistic frames and populate an
        [`EnsembleSchema`][schemas.atomistic.EnsembleSchema].

        This method initializes an [`EnsembleSchema`][schemas.atomistic.EnsembleSchema],
        calls the [`checks`][digesters.digester.Digester.checks] method to
        perform preliminary checks, prepares the inputs for digestion by calling the
        [`prepare_inputs_digester`][digesters.digester.Digester.prepare_inputs_digester]
        method, and iterates through each step of the digestion process to append
        frames to the [`EnsembleSchema`][schemas.atomistic.EnsembleSchema].

        This method keeps the entire
        [`EnsembleSchema`][schemas.atomistic.EnsembleSchema] in memory. If you have
        a large amount of data, consider using
        [`digest_chunks`][digesters.digester.Digester.digest_chunks].

        Args:
            ensemble_schema: The schema to which frames will be appended.
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            The updated schema with frames extracted from a digester.
        """
        ensemble_schema = EnsembleSchema()
        cls.checks()
        inputs_digester = cls.prepare_inputs_digester(*digester_args, **digester_kwargs)

        if not parallelize:
            for inputs_frame in cls.generate_inputs_frame(inputs_digester):
                mol_step = cls.digest_frame(inputs_frame)
                ensemble_schema.frames.append(mol_step)
        else:
            # TODO: Currently not working; the data seems out of order.
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for mol_step in executor.map(
                    cls.digest_frame, cls.generate_inputs_frame(inputs_digester)
                ):
                    ensemble_schema.frames.append(mol_step)

        return ensemble_schema

    @classmethod
    def digest_chunks(
        cls,
        digester_args: tuple[Any, ...] = tuple(),
        digester_kwargs: dict[str, Any] = dict(),
        chunk_size: int = 100,
        parallelize: bool = False,
        max_workers: int | None = None,
    ) -> Generator[EnsembleSchema, None, None]:
        """Same as [`digest`][digesters.digester.Digester.digest], but
        instead of returning a whole
        [`EnsembleSchema`][schemas.atomistic.EnsembleSchema], it will yield ones
        with a specified `chunk_size`.
        """
        cls.checks()
        inputs_digester = cls.prepare_inputs_digester(*digester_args, **digester_kwargs)

        if not parallelize:
            ensemble_schema = EnsembleSchema()
            count = 0
            for inputs_frame in cls.generate_inputs_frame(inputs_digester):
                if count > chunk_size:
                    ensemble_schema = EnsembleSchema()
                    count = 0

                mol_step = cls.digest_frame(inputs_frame)
                ensemble_schema.frames.append(mol_step)
                count += 1

                if count == chunk_size:
                    yield ensemble_schema
        else:
            # TODO: Currently not working; the data seems out of order.
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for mol_step in executor.map(
                    cls.digest_frame, cls.generate_inputs_frame(inputs_digester)
                ):
                    ensemble_schema.frames.append(mol_step)

    @classmethod
    @abstractmethod
    def get_inputs_frame(cls, inputs_digester: dict[str, Any]) -> dict[str, Any]:
        """Get the inputs for the next frame in the digestion process.

        Args:
            inputs_digester: A dictionary of inputs for the digestion process.

        Returns:
            A dictionary of inputs for the digestion process.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def next_frame(cls, inputs_digester: dict[str, Any]) -> dict[str, Any]:
        """Advance the digester inputs to the next frame in the data. This abstract
        method must be implemented in any child class as each data source may have a
        different way of advancing to the next frame.

        Args:
            inputs: A dictionary of inputs for the digestion process.

        Returns:
            A dictionary of inputs for the digestion process.
        """
        raise NotImplementedError

    @classmethod
    def generate_inputs_frame(
        cls, inputs_digester: dict[str, Any]
    ) -> Generator[dict[str, Any], None, None]:
        """Generate inputs for each frame starting from a specific frame.

        Args:
            inputs_digester: The initial inputs for the digestion process.

        Yields:
            A generator yielding input dictionaries for each frame.
        """
        # Return first frame before starting loop until the end of frames
        yield cls.get_inputs_frame(inputs_digester)
        while True:
            try:
                inputs_digester = cls.next_frame(inputs_digester)
                inputs_frame = cls.get_inputs_frame(inputs_digester)
                yield inputs_frame
            except StopIteration:
                break
            except Exception as e:
                logger.error(f"Error generating frame inputs: {e}")
                break

    @classmethod
    def digest_frame(
        cls, inputs_frame: dict[str, Any], mol_schema: MoleculeSchema | None = None
    ) -> MoleculeSchema:
        """Digest a single step of a simulation.

        This method calls each static method implemented in the child class. The static
        methods will process the data and then return a key to specify which Pydantic field
        to update and the value to update it with.

        Args:
            inputs_frame: The inputs for the MDAnalysis digester.

        Returns:
            An instance of MoleculeSchema populated with the digested data.
        """
        if mol_schema is None:
            mol_schema = MoleculeSchema()

        uuid_map = cls.get_uuid_map()
        for key_to_field, field in mol_schema.generate_fields(mol_schema):
            if not hasattr(field, "metadata") or len(field.metadata) == 0:
                logger.trace("This field does not contain any metadata")
                continue

            field_uuid = field.metadata[0].get("uuid")
            field_cadence = field.metadata[0].get("cadence")

            if field_cadence == "molecule":
                logger.debug(f"Working on molecule cadence of UUID {field_uuid}")
                if field_uuid not in uuid_map.keys():
                    logger.debug("Could not find method match of UUID")
                    continue
                method = getattr(cls, uuid_map[field_uuid])
                data = method(**inputs_frame)
                mol_schema.update({key_to_field: data})
            else:
                logger.trace("Cadence is not molecule; skipping.")

        return mol_schema
