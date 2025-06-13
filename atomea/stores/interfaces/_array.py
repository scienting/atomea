from typing import Any

from atomea.containers import AtomeaContainer
from atomea.store.interfaces import Interface
from atomea.stores import StoreKind


class ArrayInterface(Interface):
    """Provides an interface for array data."""

    def __init__(self, parent_chain: tuple[AtomeaContainer], field_name: str):
        self.parent_chain = parent_chain
        self.field_name = field_name
        self.store_kind = StoreKind.ARRAY

    def __getitem__(self, key: Any) -> Any:
        """Implement slicing with [] notation."""
        if isinstance(key, tuple):
            slices = key
        else:
            slices = (key,)

        data = self.store.read(self._path, slices=slices)
        return data

    @property
    def value(self) -> Any:
        """Get the full array without slicing."""
        data = self.store.read(self._path, slices=None)
        return data if data is not None else None
