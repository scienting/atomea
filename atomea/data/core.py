from typing import Generic

import numpy as np
import polars as pl

from atomea.data import Metadata, T, ValueOrSlice
from atomea.data.accessors import ArrayAccessor, TableAccessor
from atomea.stores import StoreKind


class Data(Generic[T]):
    """
    Descriptor that automates get/set against the appropriate store.
    Uses parent chain to access stores without storing references.
    """

    def __init__(
        self,
        dtype: type[T],
        *,
        meta: Metadata,
        default: T | None = None,
    ) -> None:
        self.dtype = dtype
        self.meta = meta
        self.default = default
        self.name: str | None = None

    def __set_name__(self, owner, name: str) -> None:
        self.name = name

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self

        # Build parent chain to root project
        parent_chain = self._get_parent_chain(obj)

        if self.meta.store is StoreKind.ARRAY:
            return ArrayAccessor(parent_chain, self.name)
        else:
            return TableAccessor(parent_chain, self.name)

    def _get_parent_chain(self, obj) -> list:
        """Build chain from current object to root project."""
        chain = []
        current = obj

        while current is not None:
            chain.append(current)
            if hasattr(current, "_stores"):  # This is the project
                break
            # Move up the hierarchy
            if hasattr(current, "_parent"):
                current = current._parent
            else:
                raise RuntimeError(f"Cannot find project root from {obj}")

        return list(reversed(chain))  # [project, ensemble, component]

    def __set__(self, obj, value: ValueOrSlice) -> None:
        if obj is None:
            raise AttributeError("Cannot set attribute on class")

        parent_chain = self._get_parent_chain(obj)
        project = parent_chain[0]

        # Determine the path based on the object hierarchy
        if len(parent_chain) == 1:
            # This is a project-level field
            path = self.name
            obj_id = getattr(obj, "id", None)
        else:
            # This is an ensemble or component-level field
            ensemble = parent_chain[1]
            path = f"{ensemble.id}/{self.name}"
            obj_id = ensemble.id

        store = project._stores[self.meta.store]

        if self.meta.store is StoreKind.ARRAY:
            if isinstance(value, tuple) and len(value) == 2:
                arr, slices = value
            else:
                arr, slices = value, None
            np_arr = np.asarray(arr, dtype=self.dtype)
            store.write(path, np_arr, slices=slices)
        else:
            # For table data
            row = {"ensemble_id": obj_id, self.name: value}
            df = pl.DataFrame([row])
            store.write(self.name, df, append=True)
