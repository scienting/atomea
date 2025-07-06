from typing import Any, Generic, Iterator

import numpy as np
import polars as pl
from loguru import logger

from atomea.containers import Container
from atomea.data import Cadence, OptionalSliceSpec, T
from atomea.stores import Store, StoreKind


class Data(Generic[T]):
    """
    Connects the user-facing containers to data in stores. Is generically typed so that
    it can interface with any store-supported data type (e.g., NumPy arrays and
    polars DataFrames).
    """

    def __init__(
        self,
        *,
        store_kind: StoreKind,
        uuid: str = "",
        description: str = "",
    ) -> None:
        self.store_kind = store_kind
        self.uuid = uuid
        self.description = description
        self.label: str | None = None
        self._container: Container | None = None
        self._parent_chain: tuple[Container, ...] = tuple()

    def __set_name__(self, owner: type, name: str) -> None:
        self.label = name

    def bind_to_container(self, container: Container) -> object:
        """Explicitly bind this Data object to a container."""
        self._container = container
        self._set_parent_chain()
        return self

    def __repr__(self):
        return f"<Data {self.label}>"

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
                logger.trace("Found root Project of {}", current)
                break
            # Move up the hierarchy
            if hasattr(current, "_parent"):
                logger.trace("Found parent container of {}", current)
                current = current._parent  # type: ignore
            else:
                raise RuntimeError(f"Cannot find project root from {self}")
        _parent_chain = tuple(reversed(chain))
        logger.trace("Setting parent chain of {} to {}", self, _parent_chain)
        self._parent_chain = _parent_chain

    @property
    def parent_chain(self) -> tuple[Container, ...]:
        """
        Chain of object references to get from the Project and any other
        `Containers` to this `Data` object. This does not include the `Data`
        object itself or the user-specified `run_id`.

        An atomea `Project` is the root container for any and all data for a
        specific project. We often nest `Containers` to group
        related information together and provide an intuitive interface. However,
        the root `Project` is the only container that keeps track of the
        `Store` backends for arrays and tables in its `Project._stores` attribute.
        This ensures that there is only one storage interface per project. All other
        nested containers contain a parent reference in `Container._parent`.

        Thus, `Data.parent_chain` traverses these `Container._parent` attributes
        backwards until it reaches the root `Project` to get access to
        `Project._stores`. However, we store it in order of `Project` to `Data` as you
        would to access this `Data` object. For example,
        `(Project, Ensemble, Microstates)`.
        """
        if len(self._parent_chain) == 0:
            self._set_parent_chain()
        return self._parent_chain

    def _get_store(self, parent_chain: tuple[Container, ...]) -> Store:
        project = parent_chain[0]
        store = project._stores[self.store_kind]  # type: ignore
        return store

    def _get_path(
        self, parent_chain: tuple[Container, ...], run_id: str | None = None
    ) -> str:
        """
        Determine the path of data in store based on the object hierarchy.

        Args:
            parent_chain:
            run_id: Optional run_id to store multiple runs for an Ensemble. This
                only impact arrays since tables store `run_id` in columns.
        """
        path = "/".join([obj.label for obj in parent_chain[1:]])
        if self.store_kind == StoreKind.ARRAY:
            if run_id is not None:
                # Cadences of ensembles should never change between run_ids.
                if self.parent_chain[-1].cadence == Cadence.MICROSTATE:
                    path += f"/{run_id}"
            path += f"/{self.label}"
        return path

    def get_store_info(self, run_id: str | None = None) -> tuple[Store, str]:
        """Determines information needed to access this data from a store.

        Returns:
            Table or Array store from the project owning this data.

            Path needed to get this data out of the store.
        """
        parent_chain = self.parent_chain
        store = self._get_store(parent_chain)
        path = self._get_path(parent_chain, run_id=run_id)
        logger.debug("Getting <Data '{}'> from <Store '{}'>", path, store.path)
        return store, path

    def write(
        self,
        data: T,
        view: OptionalSliceSpec = None,
        run_id: str | None = None,
        **kwargs: Any,
    ) -> None:
        store, path = self.get_store_info(run_id=run_id)
        store.write(path, data, view=view, **kwargs)

    def append(self, data: T, run_id: str | None = None, **kwargs: Any) -> None:
        store, path = self.get_store_info(run_id=run_id)
        store.append(path, data, **kwargs)

    def read(
        self, view: OptionalSliceSpec = None, run_id: str | None = None, **kwargs: Any
    ) -> T | None:
        store, path = self.get_store_info(run_id=run_id)
        data = store.read(path, view=view, **kwargs)
        return data  # type: ignore

    def get(
        self,
        run_id: str | None = None,
        **kwargs: Any,
    ) -> Any:
        """Get the store-specific object that represents the data stored here.

        For example, if the data is stored on disk using some type of memory map, this
        would return the memory map object, not the in-memory data. If you want to
        guarantee the data is loaded into memory, use `values`.
        """
        store, path = self.get_store_info(run_id=run_id)
        return store.get(path, **kwargs)  # type: ignore

    def iter(
        self,
        run_id: str | None = None,
        view: OptionalSliceSpec = None,
        chunk_size: int = 1,
        **kwargs: Any,
    ) -> Iterator[Any]:
        """Yield chunks of data instead of reading all into memory.

        Args:
            run_id:
            view: Slice spec for all but the first dimension.
            chunk: Number of data points of the first axis to yield at each time.
        """
        store, path = self.get_store_info(run_id=run_id)
        if store.kind != StoreKind.ARRAY:
            raise TypeError(f"Can only iterator over Array store, not {store.kind}")
        for chunk in store.iter(path, view, chunk_size, **kwargs):  # type: ignore
            yield chunk

    def next_microstate_id(self, ens_id: str, run_id: str) -> int:
        """Determine the next `microstate_id` by adding one to the currently
        largest one for this `ens_id` and `run_id`.

        Args:
            path: Path to table.
            ens_id: ID of the ensemble to query.
            run_id: ID of the run to query.

        Returns:
            The next `microstate_id`.
        """
        store, path = self.get_store_info()
        assert store.kind == StoreKind.TABLE, (
            "Can only check for microstate IDs on table data"
        )
        # Need to determine this by using the key
        df = store.query(path, ens_id=ens_id, run_id=run_id)  # type: ignore
        if df.shape == (0, 0):
            return 0

        # Get the largest microstate_id
        microstate_ids = df.get_column("microstate_id")
        microstate_id_max = microstate_ids.top_k(1).to_numpy()[0]
        return int(microstate_id_max + 1)

    def prep_dataframe(
        self, ens_id: str, run_id: str, microstate_id_next: int, data: Any
    ):
        n_microstates = data.shape[0]
        microstate_ids = np.arange(
            microstate_id_next, microstate_id_next + n_microstates
        )
        df = pl.DataFrame(
            {
                "ens_id": np.full(microstate_ids.shape, ens_id),
                "run_id": np.full(microstate_ids.shape, run_id),
                "microstate_id": microstate_ids,
                self.label: data,
            }
        )
        return df
