from typing import Any, Generic, overload

import numpy as np
import polars as pl

import atomea.typing as adt
from atomea.containers import AtomeaContainer
from atomea.data import Metadata
from atomea.data.accessors import ArrayAccessor, TableAccessor
from atomea.stores import StoreKind


class Data(Generic[adt.T]):
    """
    Data descriptor that returns exactly type T when accessed.

    The generic parameter T represents the FINAL data type that users receive,
    regardless of backend storage format.
    """

    def __init__(
        self,
        *,
        meta: Metadata,
        default: adt.T | None = None,
    ) -> None:
        self.meta = meta
        self.default = default
        self.name: str | None = None

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name

    @overload
    def __get__(self, obj: None, objtype: type) -> "Data[adt.T]": ...

    @overload
    def __get__(self, obj: Any, objtype: type | None = None) -> adt.T | None: ...

    def __get__(
        self, obj: Any | None, objtype: type | None = None
    ) -> adt.T | None | "Data[adt.T]":
        # The actual implementation of __get__
        if obj is None:
            # When accessed on the class (e.g., MyClass.descriptor_name)
            return self
        # Create accessor internally and immediately return the data
        parent_chain = self._get_parent_chain(obj)

        assert isinstance(self.name, str)

        if self.meta.store is StoreKind.ARRAY:
            accessor = ArrayAccessor[adt.T](parent_chain, self.name)
            data = accessor.value
            return data
        else:
            accessor_table = TableAccessor[adt.T](parent_chain, self.name)
            data = accessor_table.value
            return data

    def _get_parent_chain(self, obj: object) -> list[AtomeaContainer]:
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

    def __set__(self, obj: object, value: adt.ValueOrSlice[adt.T]) -> None:
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
            np_arr = np.asarray(arr)
            store.write(path, np_arr, slices=slices)
        else:
            # For table data
            row = {"ensemble_id": obj_id, self.name: value}
            df = pl.DataFrame([row])
            store.write(self.name, df, append=True)
