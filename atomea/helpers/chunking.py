from typing import Any, Iterator, Iterable


def chunker(
    chunk_size: int, n_data: int, elements: Iterable[int] | int | None = None
) -> Iterator[Any]:
    """
    Provides slice for chunking data from the first dimension.

    Args:
        chunk_size: The number of data to include in each chunk. If the total number
            of data is not evenly divisible by `chunk_size`, the final chunk will
            contain the remaining data.
        n_data: Total number of data in the first dimension.
        elements: Specific elements to chunk and yield. If `None`, then all are
            elements are chunked.

    Yields:
        Each yielded value is a subset of `data` containing up to
            `chunk_size` elements.
    """
    if elements is None:
        chunk_elements = range(0, n_data, chunk_size)
    else:
        if isinstance(elements, int):
            elements = [int]
    for start in range(0, n_data, chunk_size):
        stop = min(start + chunk_size, n_data)
        yield slice(start, stop, 1)
