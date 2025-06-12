from typing import Generic, TypeVar

import numpy as np
import numpy.typing as npt
import polars as pl

from atomea.schemas.field import FieldMeta, StoreKind

T = TypeVar("T")

SliceSpec = tuple[slice, ...] | dict[int, tuple[slice, ...]]
ValueOrSlice = T | tuple[T, SliceSpec]


class SchemaField(Generic[T]):
    """
    Descriptor that automates get/set against the appropriate store.
    """

    def __init__(
        self,
        dtype: type[T],
        *,
        meta: FieldMeta,
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

        store = obj._stores[self.meta.store]
        field_name = self.name
        obj_id = obj.id

        if self.meta.store is StoreKind.ARRAY:
            # return a getter for array data
            def getter(
                *,
                slices: tuple[slice, ...] | dict[int, tuple[slice, ...]] | None = None,
            ) -> npt.NDArray[np.generic] | None:
                path = f"{obj_id}/{field_name}"
                return store.read(path, slices=slices)

            return getter

        else:
            # return a getter for table data
            def getter(
                *,
                ensemble_id: str | None = None,
                microstate_id: int | None = None,
                filter_expr: str | None = None,
            ) -> pl.DataFrame:
                # assume the table name matches field_name
                df = store.query(
                    name=field_name,
                    ensemble_id=ensemble_id,
                    microstate_id=microstate_id,
                    filter_expr=filter_expr,
                )
                return df

            return getter

    def __set__(self, obj, value: ValueOrSlice) -> None:
        if obj is None:
            # should never happen, but satisfy mypy
            raise AttributeError("Cannot set attribute on class")

        store = obj._stores[self.meta.store]
        path = f"{obj.id}/{self.name}"

        if self.meta.store is StoreKind.ARRAY:
            if isinstance(value, tuple) and len(value) == 2:
                arr, slices = value
            else:
                arr, slices = value, None
            np_arr = np.asarray(arr, dtype=self.dtype)
            store.write(path, np_arr, slices=slices)

        else:
            # we only support full‐row writes here
            # value should be a scalar or sequence of scalars
            row = {"ensemble_id": obj.id, self.name: value}
            df = pl.DataFrame([row])
            store.write(self.name, df, append=True)
