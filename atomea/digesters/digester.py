from typing import Any, Generator

import inspect
from abc import ABC, abstractmethod
from collections.abc import Collection

from ..schemas.atomistic import EnsembleSchema, MoleculeSchema


class Digester(ABC):
    """The Digester class is an abstract base class designed to assist with digesting
    data from various computational chemistry and biology resources. This class
    provides a framework for processing, parsing, and validating data extracted from
    simulations, geometry optimizations, and other computational methods.
    """

    @classmethod
    def checks(cls):
        """Perform basic checks to raise warnings or errors before digesting.

        This method should be overridden in the child class to include specific
        checks required for the particular type of data being processed.
        """
        pass

    @classmethod
    @abstractmethod
    def prepare_step_inputs(
        cls, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        """Prepare and return the inputs necessary for the digestion process.

        This abstract method must be implemented in any child class. It should
        return a dictionary of inputs that will be used by the `digest_step` method.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            dict[str, Any]: A dictionary of inputs for the digestion process.
        """
        raise NotImplementedError

    def digest(
        self, ensemble_schema: EnsembleSchema, *args: Any, **kwargs: Collection[Any]
    ) -> EnsembleSchema:
        """Start processing, parsing, and validating data extracted from the provided simulations.

        This method calls the `checks` method to perform preliminary checks, prepares the
        inputs for digestion by calling the `prepare_step_inputs` method, and iterates through
        each step of the digestion process to append frames to the `ensemble_schema`.

        Args:
            ensemble_schema (EnsembleSchema): The schema to which frames will be appended.
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            EnsembleSchema: The updated schema with all frames extracted from a simulation.
        """
        self.checks()
        inputs = self.prepare_step_inputs(*args, **kwargs)
        for mol_step in self.digest_step(**inputs):
            ensemble_schema.frames.append(mol_step)
        return ensemble_schema

    def digest_step(
        self, *args: Any, **kwargs: Collection[Any]
    ) -> Generator[MoleculeSchema, None, None]:
        """Digest a single step of a simulation.

        This method calls each static method implemented in the child class. The static
        methods will process the data and then return a key to specify which Pydantic field
        to update and the value to update it with.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Yields:
            Generator[MoleculeSchema, None, None]: A generator yielding `MoleculeSchema` instances.
        """
        mol_schema: MoleculeSchema = MoleculeSchema()
        for name, method in inspect.getmembers(self, predicate=inspect.isfunction):
            if name not in ["digest_step", "prepare_step_inputs", "digest", "checks"]:
                _data = method(*args, **kwargs)
                mol_schema.update(_data)
        yield mol_schema
