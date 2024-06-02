from typing import Any

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from loguru import logger

from ..schemas.atomistic import EnsembleSchema


class StorageManager(ABC):
    @classmethod
    @abstractmethod
    def get_store(cls, path_file: str) -> Any:
        """Check to see if root file exists; if not, then this will initialize an
        empty file.

        Args:
            path_file: Path to file.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def initialize_group(cls, store: Any, key_group: str) -> Any:
        """Check to see if a group at `key_group` exists in `path_file`. If not,
        then it creates a group.

        Args:
            store: A store used directly from the data backend.
            key_group: Key to the group starting from root at `path_file`.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def append_array(cls, data: npt.NDArray[Any], storage: Any, key_array: str) -> None:
        """Append `data` to `key_array` located within `path_file`.

        Args:
            data: Data to append to file.
            storage: Opened file to write to.
            key_array: Key to the group starting in `storage`.
        """
        raise NotImplementedError

    @classmethod
    def _process_nested_dict(cls, nested_dict, store, prefix):
        """Recursively create groups and arrays in `store` from `nested_dict`."""
        for key, value in nested_dict.items():
            new_key = f"{prefix}/{key}"
            if isinstance(value, dict):
                # If value is a dictionary, create a new group and
                # recursively process it
                cls.initialize_group(store, new_key)
                cls._process_nested_dict(value, store, new_key)
            else:
                # Convert value to a numpy array and append it
                if isinstance(value, list):
                    value = np.array(
                        value,
                        dtype=(
                            object if isinstance(value[0], (str, dict, list)) else None
                        ),
                    )
                if value.ndim == 2:
                    value = value[None, ...]
                value = np.array(value)
                logger.debug(f"Value is: {repr(value)}")
                cls.append_array(value, store, new_key)

    @classmethod
    def process_ensemble(
        cls, ensemble_data: EnsembleSchema, store: Any, prefix: str = ""
    ) -> Any:
        """Add multiple molecules to `store` from `EnsembleSchema.frames`."""
        logger.debug("Processing ensemble")
        ensemble_dict = ensemble_data.model_dump(exclude_none=True)
        cls._process_nested_dict(ensemble_dict, store, prefix)
        return store

    @classmethod
    def digest_and_store(cls, digester, digest_kwargs, path_file):
        """Driver for digesting and writing data to storage."""
        store = cls.get_store(path_file)
        ensemble_data = digester.digest(**digest_kwargs)
        store = cls.process_ensemble(ensemble_data, store)
        return store
