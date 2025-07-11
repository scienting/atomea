# scan_engine.py
from typing import Generic, Iterator, TypeVar

from abc import ABC, abstractmethod
from enum import Enum

S = TypeVar("S", bound=Enum)


class ScanEngine(ABC, Generic[S]):
    """Abstract base class for scanning engines.

    This ABC defines the interface for any scanning engine that can
    identify patterns within a byte buffer and associate them with
    specific states. Implementations should provide methods for compiling
    patterns and streaming matches.

    Type variables:
        S: An Enum type representing the possible states a file can be in.
    """

    @abstractmethod
    def compile(self, patterns: dict[S, list[bytes]]) -> None:
        """Compile regex patterns to speed up searching.

        This method prepares the scanning engine to efficiently search for
        the provided patterns. The exact compilation mechanism will depend
        on the concrete implementation of the `ScanEngine`.

        Args:
            patterns: A dictionary where keys are state Enum members and
                values are lists of byte patterns (e.g., regex strings)
                that the engine should compile and look for.
        """

    @abstractmethod
    def stream(self, buf: bytes) -> Iterator[tuple[int, int, bytes, S]]:
        """Stream (match_start, match_end, pattern, target_state) triples.

        This method should iterate through the provided byte buffer and
        yield tuples for each match found. Each tuple should contain:

            - The inclusive starting byte offset of the match.
            - The exclusive ending byte offset of the match.
            - The actual byte pattern that was matched.
            - The target state associated with the matched pattern.

        Args:
            buf: The byte buffer to be scanned.

        Yields:
            A tuple containing the start position, end position, matched
            pattern, and the target state for each detected pattern.
        """
