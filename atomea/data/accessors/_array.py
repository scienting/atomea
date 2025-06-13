from typing import Any, Generic

import atomea.typing as adt
from atomea.data.accessors import Accessor
from atomea.stores import StoreKind


class ArrayAccessor(Accessor, Generic[adt.T]):
    """Provides property-like access and slicing for array data with proper typing."""

    @property
    def _store(self):
        """Get store from the root project."""

        project = self.parent_chain[0]
        return project._stores[StoreKind.ARRAY]  # type: ignore

    def __getitem__(self, key: Any) -> adt.T | None:
        """Implement slicing with [] notation."""
        if isinstance(key, tuple):
            slices = key
        else:
            slices = (key,)

        raw_data = self._store.read(self._path, slices=slices)
        data = self._convert_to_target_type(raw_data)
        return data

    @property
    def value(self) -> adt.T | None:
        """Get the full array without slicing."""
        data = self._store.read(self._path, slices=None)
        return data if data is not None else None

    @property
    def shape(self) -> tuple[int, ...] | None:
        """Get the shape of the array."""
        data = self._store.read(self._path, slices=None)
        return data.shape if data is not None else None
