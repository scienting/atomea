from typing import Any, Generic

from loguru import logger

from atomea.containers import AtomeaContainer
from atomea.data import Cadence, OptionalSliceSpec, T, ValueOrSlice
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
        cadence: Cadence,
        store_kind: StoreKind,
        uuid: str = "",
        description: str = "",
    ) -> None:
        self.cadence = cadence
        self.store_kind = store_kind
        self.uuid = uuid
        self.description = description
        self.name: str | None = None
        self._container: AtomeaContainer | None = None
        self._parent_chain: tuple[AtomeaContainer, ...] = tuple()

    def __set_name__(self, owner: type, name: str) -> None:
        self.name = name

    def bind_to_container(self, container: AtomeaContainer) -> object:
        """Explicitly bind this Data object to a container."""
        self._container = container
        self._set_parent_chain()
        return self

    def __repr__(self):
        return f"Data<{self.name}>"

    def _set_parent_chain(self) -> None:
        """Build chain from current object to root project and sets the attribute
        `_parent_chain`.
        """
        logger.debug("Building parent chain of {}", self)
        chain = []
        if self._container is None:
            raise RuntimeError("Data object not bound to container")

        current = self._container

        while current is not None:
            chain.append(current)
            if hasattr(current, "_stores"):  # This is the project
                logger.debug("Found root Project of {}", current)
                break
            # Move up the hierarchy
            if hasattr(current, "_parent"):
                logger.debug("Found parent container of {}", current)
                current = current._parent
            else:
                raise RuntimeError(f"Cannot find project root from {self}")
        _parent_chain = tuple(reversed(chain))
        logger.debug("Setting parent chain of {} to {}", self, _parent_chain)
        self._parent_chain = _parent_chain

    @property
    def parent_chain(self) -> tuple[AtomeaContainer, ...]:
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
            self._set_parent_chain()
        return self._parent_chain

    def get_store_info(self) -> tuple[Any, Any, Any]:
        parent_chain = self.parent_chain
        project = parent_chain[0]

        # Determine the path based on the object hierarchy
        if len(parent_chain) < 3:
            # This is a project-level field
            path = self.name
            obj_id = getattr(self, "id", None)
        else:
            # This is an ensemble or component-level field
            ensemble = parent_chain[1]
            path = f"{ensemble.id}/{self.name}"  # type: ignore
            obj_id = ensemble.id  # type: ignore

        store = project._stores[self.store_kind]  # type: ignore
        return store, path, obj_id

    def write(self, data: T, view: OptionalSliceSpec = None, **kwargs: Any) -> None:
        store, path, _ = self.get_store_info()
        store.write(path, data, view=view, **kwargs)

    def append(self, data: T, **kwargs: Any) -> None:
        store, path, _ = self.get_store_info()
        store.append(path, data, **kwargs)

    def read(self, view: OptionalSliceSpec = None, **kwargs: Any) -> T | None:
        store, path, _ = self.get_store_info()
        data = store.read(path, view=view, **kwargs)
        return data  # type: ignore

    @property
    def value(
        self,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> T | None:
        store, path, _ = self.get_store_info()
        return store.read(path, view=view, **kwargs)  # type: ignore

    def __set__(self, obj: object, value: ValueOrSlice[T]) -> None:
        if obj is None:
            raise AttributeError("Cannot set attribute on class")

        store, path, _ = self.get_store_info()

        view: OptionalSliceSpec = None
        if isinstance(value, tuple) and len(value) == 2:
            data: T = value[0]
            view = value[1]
        else:
            data = value  # type: ignore
        store.write(path, data, view=view)
