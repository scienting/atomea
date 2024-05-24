from typing import Any, Generator

import inspect
from abc import ABC, abstractmethod
from collections.abc import Collection

from ..schemas import EnsembleSchema, MoleculeSchema


class Digester(ABC):
    """Digest results into a desired atomea.

    Child classes must implement the digest method to processes all possible
    information.
    """

    @classmethod
    def checks(cls):
        pass

    @classmethod
    @abstractmethod
    def prepare_step_inputs(
        cls, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        raise NotImplementedError

    def digest(
        self, ensemble_schema: EnsembleSchema, *args: Any, **kwargs: Collection[Any]
    ) -> EnsembleSchema:
        """Digest simulations supported by [MDAnalysis](https://www.mdanalysis.org/)."""
        self.checks()
        inputs = self.prepare_step_inputs(*args, **kwargs)
        for mol_step in self.digest_step(**inputs):
            ensemble_schema.frames.append(mol_step)
        return ensemble_schema

    def digest_step(
        self, *args: Any, **kwargs: Collection[Any]
    ) -> Generator[MoleculeSchema, None, None]:
        """Digest a single step of a simulation."""
        mol_schema: MoleculeSchema = MoleculeSchema()
        for name, method in inspect.getmembers(self, predicate=inspect.isfunction):
            if name not in ["digest_step", "prepare_step_inputs", "digest", "checks"]:
                _data = method(*args, **kwargs)
                print(_data)
                mol_schema.__dict__.update(_data)
        yield mol_schema
