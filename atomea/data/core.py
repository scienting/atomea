from typing import Generic, overload

import numpy as np
import polars as pl

from atomea.containers import AtomeaContainer
from atomea.data import Metadata, T, ValueOrSlice
from atomea.store.interfaces import ArrayInterface, Interface, TableInterface
from atomea.stores import StoreKind


class Data(Generic[T]):
    """
    Connects the user-facing containers to data in stores. Is generically typed so that
    it can interface with any store-supported data type (e.g., NumPy arrays and
    polars DataFrames).
    """

    def __init__(
        self,
        *,
        meta: Metadata,
        default: T | None = None,
    ) -> None:
        self.meta = meta
        self.default = default
        self.name: str | None = None
        self._parent_chain: tuple[AtomeaContainer] = tuple()
        self._interface: ArrayInterface | TableInterface | None = None

    def _set_parent_chain(self, obj: object) -> None:
        """Build chain from current object to root project."""
        chain = []
        current = obj

        while current is not None:
            chain.append(current)
            if hasattr(current, "_stores"):  # This is the project
                break
            # Move up the hierarchy
            if hasattr(current, "_parent"):
                current = current._parent  # type: ignore
            else:
                raise RuntimeError(f"Cannot find project root from {obj}")
        self._parent_chain = tuple(reversed(chain))

    def parent_chain(self, *args):
        """
        Chain of object references to get from the Project and any other
        `AtomeaContainers` to this `Data` object.

        An atomea `Project` is the root container for any and all data for a
        specific project. We often nest `AtomeaContainers` to group
        related information together and provide an intuitive interface. However,
        the root `Project` is the only container that keeps track of the
        `Store` backends for arrays and tables in its `Project._stores` attribute.
        This ensures that there is only one storage interface per project. All other
        nested containers contain a parent reference in `AtomeaContainer._parent`.

        Thus, `Data.parent_chain` traverses these `AtomeaContainer._parent` attributes
        backwards until it reaches the root `Project` to get access to
        `Project._stores`. However, we store it in order of `Project` to `Data` as you
        would to access this `Data` object. For example,
        `(Project, Ensemble, Microstates)`.
        """
        if len(self._parent_chain) == 0:
            self._set_parent_chain(self, *args)
        return self._parent_chain

    def interface(self, obj: object) -> ArrayInterface | TableInterface:
        """Interface for this Data. Cannot changed after it is called once."""
        if self._interface:
            return self._interface

        parent_chain = self.parent_chain(obj)

        assert isinstance(self.name, str)

        if self.meta.store is StoreKind.ARRAY:
            interface = ArrayInterface(parent_chain, self.name)
        elif self.meta.store is StoreKind.TABLE:
            interface = TableInterface(parent_chain, self.name)
        else:
            raise TypeError("Unknown store type of {}", self.meta.store)
        return interface

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name

    @overload
    def __get__(self, obj: None, objtype: type) -> "Data[T]": ...

    @overload
    def __get__(self, obj: object, objtype: type | None = None) -> T | None: ...

    def __get__(
        self, obj: object | None, objtype: type | None = None, *args, **kwargs
    ) -> T | None | "Data[T]":
        if obj is None:
            return self
        interface = self.interface(obj)
        return interface.__getitem__(*args, **kwargs)

    def __set__(self, obj: object, value: ValueOrSlice[T]) -> None:
        if obj is None:
            raise AttributeError("Cannot set attribute on class")

        parent_chain = self.parent_chain(obj)
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

        store = project.stores[self.meta.store]  # type: ignore

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
