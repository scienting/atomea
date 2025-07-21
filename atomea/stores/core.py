# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

import os
from abc import ABC, abstractmethod
from pathlib import Path

import atomea.typing as adt
from atomea.data import OptionalSliceSpec
from atomea.stores import DiskFormat, StoreKind


class Store(ABC):
    """
    An abstract interface for handling data stored in memory or on disk.
    """

    kind: StoreKind
    """Which [kind of store][stores.kind.StoreKind] this is."""

    _store: Any
    """Underlying interface for accessing data using this store."""

    def __init__(
        self,
        path: Path | str = "",
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            path: Path where data will be stored on disk. We will create this
                directory if it does not exist depending on the value of
                `disk_format`.
            mode: Persistence modes:

                - `r` means read only (must exist);
                - `r+` means read/write (must exist);
                - `a` means read/write (create if doesn't exist);
                - `w` means create (overwrite if exists);
                - `w-` means create (fail if exists).
            disk_format: File format when writing data to disk. If left as the
                default `DiskFormat.NONE`, then we will not check for the existence
                of `path`.
        """
        self.path = path
        self.mode = mode
        self.disk_format = disk_format
        if path != "" and disk_format != DiskFormat.NONE:
            os.makedirs(path, exist_ok=True)

    @abstractmethod
    def write(
        self,
        path: Path | str,
        data: Any,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        """Write or update data in the underlying store.

        This method is responsible for writing new data or updating existing data
        at the specified `path`. The exact behavior (e.g., atomic overwrite,
        partial update) depends on the `view` argument and the capabilities
        of the underlying data storage mechanism.

        For array-like data (e.g., NumPy, Zarr), `data` should be an array-like
        object. If `view` is provided, it will attempt to write the `data` into
        the specified slice of the existing array. If `view` is `None`, it will
        typically overwrite or create the array at `path`.

        For table-like data (e.g., Polars DataFrames), `data` should be a
        DataFrame or Series. If `view` is `None`, it will typically write the
        entire DataFrame to the `path`. If `view` is used to specify a subset
        (e.g., a column name or row selection), it will attempt to update that
        specific portion of the table.

        Args:
            path: The file path or identifier where the data should be written or
                updated.
            data: The data to write. This can be a NumPy array, Polars DataFrame,
                  Polars Series, or other compatible data structure depending on
                  the data type being managed.
            view: An optional slice specification for partial updates. For arrays,
                  this defines the region to write `data` into. For tables, this
                  might specify rows or columns to update. If `None`, the entire
                  data at `path` is typically overwritten or created.
            **kwargs: Additional keyword arguments specific to the underlying
                  storage mechanism (e.g., compression options, encoding,
                  write modes for Polars).

        Note:
            This method updates the internal state of the data store, but it **does not
            guarantee immediate persistence to disk**. For some backends (like Zarr),
            writes might be flushed automatically or batched. For others (like
            memory-mapped NumPy files or certain table formats), changes might remain
            in a buffer or cache until explicitly flushed. To ensure data is written
            to durable storage, you **must** call the `flush()` method.
        """

    def flush(self, **kwargs: Any) -> None:
        """Force all pending data changes to be written to disk.

        This method ensures that any modifications made via `write()` or other
        operations are explicitly synchronized with the underlying durable storage
        (e.g., hard drive). It is crucial for guaranteeing data integrity and
        persistence, especially in scenarios where immediate durability is required
        or where the underlying storage mechanism employs buffering or lazy writing.

        For certain backends, such as Zarr, writes are often atomic and immediately
        flushed to disk, making this `flush` method a no-op or a minimal operation.
        However, for other data types like memory-mapped NumPy arrays or Polars
        DataFrames written to files, calling `flush()` will trigger the actual
        write operation from memory buffers to the physical disk.

        It is recommended to call this method after a series of `write` operations
        or before exiting an application to prevent data loss.

        Args:
            **kwargs: Additional keyword arguments that might be passed to the
                underlying flushing mechanism (e.g., specific sync options).
        """

    @abstractmethod
    def append(
        self,
        path: Path | str,
        data: Any,
        **kwargs: Any,
    ) -> None:
        """Append data to the existing data structure at the given path.

        This method adds `data` to the end of the existing dataset at `path`.
        Its behavior is type-specific:

        For array-like data (e.g., NumPy, Zarr):
        The `data` provided will be concatenated along a predefined axis (e.g.,
        the first dimension, representing new records or entries) to the existing
        array at `path`. The shape of `data` must be compatible with the existing
        array for the append operation to succeed. This operation might implicitly
        resize the underlying storage.

        For table-like data (e.g., Polars DataFrames):
        The `data` (which should be a DataFrame or an iterable of records/rows)
        will be appended as new rows to the table at `path`. The columns of the
        appended `data` must match or be compatible with the existing table's schema.

        Args:
            path: The file path or identifier of the data structure to append to.
            data: The data to append. For arrays, this should be an array-like
                  object. For tables, this should typically be a DataFrame or
                  a sequence of records that can be converted into new rows.
            **kwargs: Additional keyword arguments specific to the underlying
                      storage mechanism (e.g., options for column matching
                      during table appends, or handling of missing columns).

        Notes:
            Similar to `write()`, this method updates the data store, but it **does
            not guarantee immediate persistence to disk**. Changes may be buffered
            internally. To ensure data is durably saved to disk, you **must** call
            the `flush()` method after appending.
        """

    @abstractmethod
    def get(
        self,
        path: Path | str,
        **kwargs: Any,
    ) -> Any | None:
        """Access a reference of the data at the given path.

        This method provides a way to interact with the data without necessarily
        loading its entire contents into active memory. It's ideal for large
        datasets where you might only need to access parts of the data, or when
        working with memory-mapped files or lazy-loading structures.

        For arrays, this might return a Zarr array, a NumPy memory-mapped array,
        or an in-memory NumPy array if the data is small or already loaded.
        For tables, it will return a DataFrame, which
        might internally manage memory efficiently for large datasets or provide
        a view into a file-backed store.

        Use this method when you want to avoid immediate, full memory allocation
        and prefer a more performant way to access potentially large datasets.

        Args:
            path: The file path or identifier where the data is located.
            **kwargs: Additional keyword arguments to pass to the underlying
                      data loading mechanism (e.g., specific Zarr or Pandas
                      options).

        Returns:
            A data object (e.g., Zarr array, memory-mapped NumPy array, Pandas DataFrame)
                that represents the data. This object may provide a "view" or a reference,
                and not necessarily a fully in-memory copy of all data. Returns None if
                the data cannot be read or found.
        """

    @abstractmethod
    def read(
        self, path: Path | str, view: OptionalSliceSpec = None, **kwargs: Any
    ) -> adt.OptionalPassableData:
        """Retrieve the data at the given path, guaranteed to be loaded into memory.

        This method ensures that the requested data (or a specified slice/series of it)
        is fully loaded into Python's active memory as a standard in-memory object.
        This is suitable when you need to perform operations that require all data
        to be immediately accessible in RAM, such as complex computations,
        modifications, or when the dataset size is manageable for in-memory processing.

        For arrays, this will always return a standard NumPy array with all its
        contents loaded into memory. For tables, if no `view` (slice) or specific
        `column_name` is provided via `kwargs`, this method will return the
        entire DataFrame fully loaded into memory. If a `column_name` is specified
        or implied via `view`, it will return an in-memory Series.

        Args:
            path: The file path or identifier where the data is located.
            view: An optional slice specification to retrieve a subset of the
                data. For arrays, this could be a standard NumPy slice. For
                tables, this might implicitly select rows or specific columns
                if combined with other arguments (see `**kwargs`).
            **kwargs: Additional keyword arguments. For tables, this may include
                a `column_name` (str) to retrieve a specific series, or
                `filters` (dict/list) to apply before loading into memory.
                Other arguments will be passed to the underlying loading
                mechanism.

        Returns:
            The requested data as a fully in-memory object (e.g., NumPy array,
            Pandas DataFrame, or Pandas Series). Returns None if the data
            cannot be retrieved or found.
        """

    @abstractmethod
    def available(self) -> list[str]:
        """List the names or paths of all data entries currently managed by this store.

        This method provides an overview of the data accessible through the current
        data manager instance. The returned names typically correspond to the `path`
        arguments used in `read()`, `get()`, `write()`, and `append()` methods,
        allowing users to discover what data is available.

        The exact format of the names (e.g., file paths, internal identifiers)
        depends on the specific implementation of the data manager and the
        underlying storage system.

        Returns:
            A list of strings, where each string represents the name or path
                of an available data entry within the store.
        """
