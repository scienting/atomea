from typing import Generic, cast

import atomea.typing as adt
from atomea.data.accessors import Accessor
from atomea.stores import StoreKind


class TableAccessor(Accessor, Generic[adt.T]):
    """Provides property-like access for table data with proper typing."""

    @property
    def _store(self):
        """Get store from the root project."""
        project = self.parent_chain[0]
        return project._stores[StoreKind.TABLE]  # type: ignore

    def __getitem__(
        self, key: slice | str | int | tuple[str | int | None, ...]
    ) -> adt.T:
        """Support bracket notation with positional arguments.

        Usage:
            accessor[:]                           # All data
            accessor["e1"]                        # Ensemble "e1"
            accessor["e1", 0]                     # Ensemble "e1", microstate 0
            accessor["e1", 0, "energy > 100"]     # All three parameters
        """
        ensemble_id = None
        microstate_id = None
        filter_expr = None

        if key == slice(None):
            # accessor[:] - get all data
            pass
        elif isinstance(key, (str, int)):
            # accessor["e1"] or accessor[0] - single argument
            if isinstance(key, str):
                ensemble_id = key
            else:
                microstate_id = key
        elif isinstance(key, tuple):
            # accessor["e1", 0] or accessor["e1", 0, "filter"] - multiple arguments
            if len(key) >= 1 and key[0] is not None:
                ensemble_id = str(key[0])
            if len(key) >= 2 and key[1] is not None:
                microstate_id = int(key[1])
            if len(key) >= 3 and key[2] is not None:
                filter_expr = str(key[2])

        df = self._store.query(
            name=self.field_name,
            ensemble_id=ensemble_id,
            microstate_id=microstate_id,
            filter_expr=filter_expr,
        )
        return cast(adt.T, df)

    @property
    def value(self) -> adt.T:
        """Get the full table data."""
        df = self._store.query(name=self.field_name)
        return cast(adt.T, df)
