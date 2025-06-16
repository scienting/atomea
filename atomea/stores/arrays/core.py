from typing import Any

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt

from atomea.data import OptionalSliceSpec
from atomea.stores import ArrayDiskFormats, DiskFormat, Store


class ArrayStore(Store, ABC):
    """
    Abstract interface for storing and retrieving per-ensemble and per-microstate arrays,
    e.g., coordinates, velocities, forces, etc.
    """

    def __init__(
        self,
        path: str,
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            path: Path to directory where arrays will be stored. For example, this
                should be to a Zarr store ending in `.zarr` or a directory of
                `<prj_name>.arrays` for NumPy.
            disk_format: File format when writing arrays to disk.
        """
        assert disk_format in ArrayDiskFormats
        super().__init__(path, disk_format=disk_format, **kwargs)

    @abstractmethod
    def create(
        self,
        path: str,
        shape: tuple[int, ...],
        overwrite: bool = False,
        **kwargs: Any,
    ) -> Any: ...

    @abstractmethod
    def write(
        self,
        path: str,
        data: Any,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> None:
        """
        Write a new array at the given path. Overwrites if it exists.

        Args:
            path: logical path or key within the ensemble (e.g., "coords/0").
            array: numpy array to store.
        """
        ...

    @abstractmethod
    def append(self, path: str, data: npt.NDArray[np.generic], **kwargs: Any) -> None:
        """
        Append data to an existing array at the given path, creating it if needed.
        This is useful for adding new microstates one by one.

        Args:
            path: logical path or key within the ensemble.
            array: numpy array to append along the first axis.
        """
        ...

    @abstractmethod
    def read(
        self,
        path: str,
        view: OptionalSliceSpec = None,
        **kwargs: Any,
    ) -> npt.NDArray[np.generic] | None:
        """
        Retrieve the full array stored at the given path.

        Args:
            path: logical path or key within the ensemble.
            view: either a tuple of slice objects corresponding to each dimension,
                    or a dict of {axis_index: slice} for partial reads.
                    If None, returns the full array.
        Returns:
            numpy array previously stored.
        """
        ...

    @abstractmethod
    def available(self) -> list[str]:
        """
        List all array paths/keys currently stored in this ensemble.

        Returns:
            A list of logical paths (e.g., ["coords", "velocities"]).
        """
        ...
