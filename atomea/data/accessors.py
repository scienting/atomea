import numpy as np
import numpy.typing as npt
import polars as pl

from atomea.data import SliceSpec
from atomea.stores import StoreKind


class ArrayAccessor:
    """Provides property-like access and slicing for array data."""

    def __init__(self, parent_chain: list, field_name: str):
        self.parent_chain = parent_chain  # [project, ensemble, component]
        self.field_name = field_name

    @property
    def _store(self):
        """Get store from the root project."""
        project = self.parent_chain[0]
        return project._stores[StoreKind.ARRAY]

    @property
    def _path(self):
        """Build path from parent chain."""
        ensemble = self.parent_chain[1] if len(self.parent_chain) > 1 else None
        if ensemble is None:
            # This is a project-level field
            return self.field_name
        return f"{ensemble.id}/{self.field_name}"

    def __getitem__(self, key) -> npt.NDArray[np.generic] | None:
        """Support slicing with [] notation."""
        if isinstance(key, tuple):
            slices = key
        else:
            slices = (key,)
        return self._store.read(self._path, slices=slices)

    def __call__(
        self, *, slices: SliceSpec | None = None
    ) -> npt.NDArray[np.generic] | None:
        """Support the original function call interface."""
        return self._store.read(self._path, slices=slices)

    @property
    def value(self) -> npt.NDArray[np.generic] | None:
        """Get the full array without slicing."""
        return self._store.read(self._path, slices=None)

    @property
    def shape(self) -> tuple[int, ...] | None:
        """Get the shape of the array."""
        data = self._store.read(self._path, slices=None)
        return data.shape if data is not None else None

    @property
    def dtype(self) -> np.dtype | None:
        """Get the dtype of the array."""
        data = self._store.read(self._path, slices=None)
        return data.dtype if data is not None else None


class TableAccessor:
    """Provides property-like access for table data."""

    def __init__(self, parent_chain: list, field_name: str):
        self.parent_chain = parent_chain
        self.field_name = field_name

    @property
    def _store(self):
        """Get store from the root project."""
        project = self.parent_chain[0]
        return project._stores[StoreKind.TABLE]

    def __call__(
        self,
        *,
        ensemble_id: str | None = None,
        microstate_id: int | None = None,
        filter_expr: str | None = None,
    ) -> pl.DataFrame:
        """Support the original function call interface."""
        return self._store.query(
            name=self.field_name,
            ensemble_id=ensemble_id,
            microstate_id=microstate_id,
            filter_expr=filter_expr,
        )

    @property
    def value(self) -> pl.DataFrame:
        """Get the full table data."""
        return self._store.query(name=self.field_name)
