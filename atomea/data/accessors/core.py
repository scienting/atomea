from typing import Any, TypeAlias, Union, cast

from abc import ABC

import numpy as np

import atomea.typing as adt
from atomea.containers import AtomeaContainer


class Accessor(ABC):
    def __init__(self, parent_chain: list[AtomeaContainer], field_name: str):
        self.parent_chain = parent_chain
        self.field_name = field_name

    @staticmethod
    def _get_numpy_type(type_target: TypeAlias) -> np.dtype:
        """
        Get numpy data type given the target type. Defaults to `np.float64`.
        """
        # Handle the type_target which might be a union (T | None)
        if hasattr(type_target, "__origin__") and type_target.__origin__ is Union:
            # Extract the non-None type from Optional[T]
            args = type_target.__args__
            non_none_type = next((arg for arg in args if arg is not type(None)), None)
            if non_none_type:
                return adt.DTYPE_TO_NUMPY.get(non_none_type, np.float64)  # type: ignore

        return adt.DTYPE_TO_NUMPY.get(type_target, np.float64)  # type: ignore

    @property
    def _path(self):
        """Build path from parent chain."""
        ensemble: AtomeaContainer | None = (
            self.parent_chain[1] if len(self.parent_chain) > 1 else None
        )
        if ensemble is None:
            # This is a project-level field
            return self.field_name
        else:
            return f"{ensemble.id}/{self.field_name}"  # type: ignore

    def _convert_to_target_type(self, raw_data: Any) -> adt.T | None:
        """Convert raw backend data to target type T."""
        if raw_data is None:
            return None

        # Convert to numpy array with correct dtype
        try:
            np_array = np.asarray(raw_data)
            return cast(adt.T, np_array)
        except (ValueError, TypeError):
            # Handle conversion errors gracefully
            return None
