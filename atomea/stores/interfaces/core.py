from abc import ABC

from atomea.containers import AtomeaContainer
from atomea.stores import StoreKind


class Interface(ABC):
    """Interface for interfaces."""

    store_kind: StoreKind

    def __init__(self, parent_chain: tuple[AtomeaContainer], field_name: str):
        self.parent_chain = parent_chain
        self.field_name = field_name

    @property
    def store(self):
        """Get store from the root project."""

        project = self.parent_chain[0]
        return project._stores[self.store_kind]  # type: ignore

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

    @property
    def shape(self) -> tuple[int, ...] | None:
        """Get the shape of the array or dataframe."""
        data = self.store.read(self._path, slices=None)
        return data.shape if data is not None else None
