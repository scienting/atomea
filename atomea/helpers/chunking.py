# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

from collections.abc import Iterable, Iterator

from loguru import logger


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
            elements are chunked. If provided, elements are yielded in the order they
            appear in the input iterable, without de-duplication or sorting.

    Yields:
        Each yielded value is a subset of `data` containing up to
            `chunk_size` elements.
    """
    logger.debug(
        "Generating chunk slices of size {} for {} elements", chunk_size, n_data
    )
    if elements is None or elements == slice(None, None, None):
        for start in range(0, n_data, chunk_size):
            stop = min(start + chunk_size, n_data)
            yield slice(start, stop, 1)
    else:
        # Handle single integer element uniformly as an iterable
        elements_iterable: Iterable[int]
        if isinstance(elements, int):
            elements_iterable = (elements,)
        else:
            elements_iterable = elements

        current_chunk = []
        for element in elements_iterable:
            if not (0 <= element < n_data):
                raise RuntimeError(
                    f"Element {element} is out of bounds! Must be between 0 and {n_data - 1}."
                )

            current_chunk.append(element)
            if len(current_chunk) == chunk_size:
                yield current_chunk
                current_chunk = []

        # Yield any remaining elements in the last chunk
        if current_chunk:
            yield current_chunk
