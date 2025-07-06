from typing import Any, Iterator

from abc import ABC
from pathlib import Path

from atomea.data import OptionalSliceSpec
from atomea.helpers import chunker
from atomea.stores import ArrayDiskFormats, DiskFormat, Store, StoreKind


class ArrayStore(Store, ABC):
    """
    Abstract interface for storing and retrieving arrays,
    e.g., coordinates, velocities, forces.
    """

    kind = StoreKind.ARRAY

    def __init__(
        self,
        path: Path | str,
        mode: str = "r",
        disk_format: DiskFormat = DiskFormat.NONE,
        **kwargs: Any,
    ) -> None:
        assert disk_format in ArrayDiskFormats or disk_format == DiskFormat.NONE
        super().__init__(path, mode=mode, disk_format=disk_format, **kwargs)

    def iter(
        self,
        path: Path | str,
        view: OptionalSliceSpec = None,
        chunk_size: int = 1,
        **kwargs: Any,
    ) -> Iterator[Any]:
        """Yield chunks of data instead of reading all into memory.

        Args:
            path:
            view: View for all but the first dimension.
            chunk: Number of data points of the first axis to yield.
        """
        z = self.get(path)
        if z is None:
            return None
        n_items = z.shape[0]
        for chunk in chunker(chunk_size, n_items):
            if view is None:
                _view = (chunk,)
            else:
                _view = (chunk, *view)  # type: ignore
            yield z.get_orthogonal_selection(_view)
