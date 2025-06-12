"""Help with generic types for Data"""

from typing import TypeVar

T = TypeVar("T")

SliceSpec = tuple[slice, ...] | dict[int, tuple[slice, ...]]
ValueOrSlice = T | tuple[T, SliceSpec]
