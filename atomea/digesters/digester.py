from typing import Any, Collection, Generator, Literal

import inspect
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor

from loguru import logger

from ..schemas.atomistic import EnsembleSchema


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
        return None

    @classmethod
    def get_uuid_map(cls):
        """
        Update the function UUID map by inspecting the class methods
        decorated with [`@SchemaUUID`][digesters.ids.SchemaUUID].

        This method scans through all the methods in the class, identifies those
        decorated with [`@SchemaUUID`][digesters.ids.SchemaUUID], and constructs a
        dictionary mapping the UUIDs to the method names. This map is used to
        dynamically call methods based on their UUIDs during the data digestion process.

        By using [`@SchemaUUID`][digesters.ids.SchemaUUID], each method that processes
        a part of the input data can be easily identified and called based on its UUID.
        This allows for a flexible and dynamic way to handle various data processing
        tasks, ensuring that each piece of data is processed by the appropriate method.

        Returns:
            A dictionary mapping UUIDs to method names.

        Raises:
            `RuntimeError` if there are duplicate UUIDs.

        Example:
            If a method called `coordinates` is decorated with
            @SchemaUUID("81c7cec9-beec-4126-b6d8-91bee28951d6"), the returned
            dictionary will include an entry:
            `{"81c7cec9-beec-4126-b6d8-91bee28951d6": "coordinates"}`

        Notes:
            This method only includes methods that

            -   are callable,
            -   do not start with `__`, and
            -   have the `__uuid__` attribute.
        """
        uuid_map: dict[str, str] = {}
        for name, method in inspect.getmembers(cls, predicate=inspect.isfunction):
            if name[:2] == "__":
                continue
            if callable(method) and hasattr(method, "__uuid__"):
                uuid = method.__uuid__
                if uuid in uuid_map.keys():
                    raise RuntimeError(f"Found duplicate UUID: {uuid}")
                uuid_map[method.__uuid__] = name
        return uuid_map

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
        ensemble_schema: EnsembleSchema | None = None,
        digester_args: tuple[Any, ...] | None = None,
        digester_kwargs: dict[str, Any] | None = None,
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
            digester_args: Arguments to pass into
                [`prepare_inputs_digester`][digesters.digester.Digester.prepare_inputs_digester].
            digester_kwargs: Keyword arguments to pass into
                [`prepare_inputs_digester`][digesters.digester.Digester.prepare_inputs_digester].
            parallelize: Execute concurrently.
            max_workers: Maximum number of workers for concurrent operation.

        Returns:
            The updated schema with frames extracted from a digester.
        """
        if ensemble_schema is None:
            ensemble_schema = EnsembleSchema()
        if digester_args is None:
            digester_args = tuple()
        if digester_kwargs is None:
            digester_kwargs = {}

        cls.checks()
        inputs_digester = cls.prepare_inputs_digester(*digester_args, **digester_kwargs)

        schema_map = ensemble_schema.get_schema_map(ensemble_schema)
        mol_index = 0
        # Digest all frames with a cadence of molecule
        if not parallelize:
            for inputs_frame in cls.gen_inputs_frame(inputs_digester):
                mol_data = cls.digest_frame(inputs_frame, schema_map)
                mol_index = ensemble_schema.update_atomistic(
                    mol_data, schema_map=schema_map, mol_index=mol_index
                )
        else:
            logger.error("Parallel operation is not yet supported.")
            # TODO: Currently not working; the data seems out of order.
            # with ThreadPoolExecutor(max_workers=max_workers) as executor:
            #     for mol_frame in executor.map(
            #         cls.digest_frame, cls.gen_inputs_frame(inputs_digester)
            #     ):
            #         ensemble_schema.frames.append(mol_frame)

        # Cleanup molecule arrays, must be done before ensemble.
        ensemble_schema._trim_molecule_arrays(mol_index, schema_map)

        # Digest all ensemble-cadence properties using the last frame
        ensemble_data = cls.digest_frame(
            inputs_frame, schema_map, cadence_eval="ensemble"
        )
        ensemble_schema.update_atomistic(ensemble_data, schema_map=schema_map)

        return ensemble_schema

    @classmethod
    def digest_chunks(
        cls,
        ensemble_schema: EnsembleSchema | None = None,
        digester_args: tuple[Any, ...] | None = None,
        digester_kwargs: dict[str, Any] | None = None,
        chunk_size: int = 100,
        parallelize: bool = False,
        max_workers: int | None = None,
    ) -> Generator[EnsembleSchema, None, None]:
        """Same as [`digest`][digesters.digester.Digester.digest], but
        instead of returning a whole
        [`EnsembleSchema`][schemas.atomistic.EnsembleSchema], it will generate ones
        with a specified `chunk_size`.

        Args:
            ensemble_schema: The schema to which frames will be appended.
            digester_args: Arguments to pass into
                [`prepare_inputs_digester`][digesters.digester.Digester.prepare_inputs_digester].
            digester_kwargs: Keyword arguments to pass into
                [`prepare_inputs_digester`][digesters.digester.Digester.prepare_inputs_digester].
            chunk_size: Number of frames to process before yielding an
                [`EnsembleSchema`][schemas.atomistic.EnsembleSchema].
            parallelize: Execute concurrently.
            max_workers: Maximum number of workers for concurrent operation.
        """
        if ensemble_schema is None:
            ensemble_schema = EnsembleSchema()
        ensemble_schema_orig = ensemble_schema.copy()
        if digester_args is None:
            digester_args = tuple()
        if digester_kwargs is None:
            digester_kwargs = {}

        cls.checks()
        inputs_digester = cls.prepare_inputs_digester(*digester_args, **digester_kwargs)

        schema_map = ensemble_schema.get_schema_map(ensemble_schema)
        mol_index = 0
        if not parallelize:
            ensemble_schema = ensemble_schema_orig.copy()
            count = 0
            for inputs_frame in cls.gen_inputs_frame(inputs_digester):
                if count > chunk_size:
                    ensemble_schema = ensemble_schema_orig.copy()
                    count = 0

                mol_data = cls.digest_frame(inputs_frame, schema_map)
                mol_index = ensemble_schema.update_atomistic(
                    mol_data, schema_map=schema_map, mol_index=mol_index
                )
                count += 1

                # Cleanup molecule arrays, must be done before ensemble.
                ensemble_schema._trim_molecule_arrays(mol_index, schema_map)

                # Digest all ensemble-cadence properties using the last frame
                ensemble_data = cls.digest_frame(
                    inputs_frame, schema_map, cadence_eval="ensemble"
                )
                ensemble_schema.update_atomistic(ensemble_data, schema_map=schema_map)

                if count == chunk_size:
                    yield ensemble_schema
        else:
            logger.error("Parallel operation is not yet supported.")
            # TODO: Currently not working; the data seems out of order.
            # with ThreadPoolExecutor(max_workers=max_workers) as executor:
            #     for mol_frame in executor.map(
            #         cls.digest_frame, cls.gen_inputs_frame(inputs_digester)
            #     ):
            #         ensemble_schema.frames.append(mol_frame)

    @classmethod
    @abstractmethod
    def get_inputs_frame(cls, inputs_digester: dict[str, Any]) -> dict[str, Any]:
        """Builds dictionary of keyword arguments for the current frame specified
        in `inputs_digester`. This is called for every frame.

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
    def gen_inputs_frame(
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
        cls,
        inputs_frame: dict[str, Any],
        schema_map: dict[str, dict[str, str]],
        cadence_eval: Literal["molecule", "ensemble"] = "molecule",
    ) -> dict[str, Any]:
        """
        Digest a single frame of input data into a
        [`MoleculeSchema`][schemas.atomistic.MoleculeSchema].

        This method processes a single frame of data by invoking static methods
        implemented in the child digester class. These static methods with a
        [`SchemaUUID`][digesters.ids.SchemaUUID] are responsible
        for processing specific parts of the frame input data and returning
        key-value pairs that correspond to fields in the
        [`MoleculeSchema`][schemas.atomistic.MoleculeSchema].

        Args:
            inputs_frame: The inputs for the frame digestion process.
                This dictionary should contain all necessary data for processing a
                single frame.
            schema_map: A mapping of UUIDs to field keys from
                [`get_schema_map`][schemas.atomistic.ensemble.EnsembleSchema.get_schema_map]
            cadence_eval: Cadence of properties to evaluate and digest.

        Returns:
            Data parsed or computed for this frame. Keys are field keys and values
            are the data from this frame.

        Raises:
            AttributeError: If the static method corresponding to a field's UUID
                is not found in the class.
            Exception: For any other exceptions that occur during the processing
                of the frame.

        Notes:
            - The method relies on metadata defined within the fields of MoleculeSchema
              to determine which static method to call for processing each field.
            - Each field in the MoleculeSchema should have metadata that includes a
              'uuid' and optionally a 'cadence'. The 'cadence' should be set to 'molecule'
              to indicate that the field is processed per molecule.
            - Static methods in the child class should be decorated with @SchemaUUID
              to associate them with the corresponding fields in MoleculeSchema.

        Example:
            Suppose `inputs_frame` contains data for atomic coordinates, the static
            method decorated with the appropriate UUID will be called to process these
            coordinates, and the resulting values will be assigned to the corresponding
            field in MoleculeSchema.

        """
        digester_map = cls.get_uuid_map()
        mol_data = {}
        for field_uuid, field_info in schema_map.items():
            field_cadence = field_info["cadence"]
            field_key = field_info["field_key"]

            if field_cadence == cadence_eval:
                logger.debug(f"Working on {cadence_eval} cadence of UUID {field_uuid}")
                if field_uuid not in digester_map.keys():
                    logger.debug("Could not find method match of UUID")
                    continue
                method = getattr(cls, digester_map[field_uuid])
                data = method(**inputs_frame)
                mol_data[field_key] = data
            else:
                logger.trace(f"Cadence is not {cadence_eval}; skipping.")

        return mol_data
