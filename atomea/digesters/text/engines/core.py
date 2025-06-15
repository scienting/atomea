# scan_engine.py
from typing import Generic, Iterator, TypeVar

from abc import ABC, abstractmethod
from enum import Enum

S = TypeVar("S", bound=Enum)


class ScanEngine(ABC, Generic[S]):
    """Return (match_start, match_end, target_state) triples."""

    @abstractmethod
    def compile(self, patterns: dict[S, list[bytes]]) -> None: ...

    @abstractmethod
    def stream(self, buf: bytes) -> Iterator[tuple[int, int, str, S]]: ...
