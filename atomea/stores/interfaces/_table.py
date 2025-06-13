from typing import Any

from atomea.containers import AtomeaContainer
from atomea.store.interfaces import Interface
from atomea.stores import StoreKind


class TableInterface(Interface):
    """Provides an interface for table data."""

    def __init__(self, parent_chain: tuple[AtomeaContainer], field_name: str):
        self.parent_chain = parent_chain
        self.field_name = field_name
        self.store_kind = StoreKind.TABLE

    def __getitem__(self, key: slice | str | int | tuple[str | int | None, ...]) -> Any:
        """Support bracket notation with positional arguments.

        Usage:
            interface[:]                           # All data
            interface["e1"]                        # Ensemble "e1"
            interface["e1", 0]                     # Ensemble "e1", microstate 0
            interface["e1", 0, "energy > 100"]     # All three parameters
        """
        ensemble_id = None
        microstate_id = None
        filter_expr = None

        if key == slice(None):
            # interface[:] - get all data
            pass
        elif isinstance(key, (str, int)):
            # interface["e1"] or interface[0] - single argument
            if isinstance(key, str):
                ensemble_id = key
            else:
                microstate_id = key
        elif isinstance(key, tuple):
            # interface["e1", 0] or interface["e1", 0, "filter"] - multiple arguments
            if len(key) >= 1 and key[0] is not None:
                ensemble_id = str(key[0])
            if len(key) >= 2 and key[1] is not None:
                microstate_id = int(key[1])
            if len(key) >= 3 and key[2] is not None:
                filter_expr = str(key[2])

        df = self.store.query(
            name=self.field_name,
            ensemble_id=ensemble_id,
            microstate_id=microstate_id,
            filter_expr=filter_expr,
        )
        return df

    @property
    def value(self) -> Any:
        """Get the full table data."""
        df = self.store.query(name=self.field_name)
        return df
