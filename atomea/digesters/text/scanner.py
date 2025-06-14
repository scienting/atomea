from typing import Generic, TypeVar

from abc import ABC, abstractmethod
from dataclasses import dataclass

from atomea.digesters.text import Pattern, StateTransition, TransitionHint

T = TypeVar("T")


@dataclass
class ScannerConfig:
    chunk_size: int = 65536
    """Size of chunk in bytes. Defaults to 64KB"""
    overlap_size: int = 1024
    """Size of overlap between chunks in bytes. Defaults to 1KB"""
    max_lookahead: int = 4096
    parallel_threshold: int = 10_000_000  # 10MB


class StateScanner(ABC, Generic[T]):
    """Abstract scanner for identifying state transitions"""

    @abstractmethod
    def get_patterns(self) -> dict[T, list[Pattern]]:
        """Return patterns that trigger transitions TO each state"""
        pass

    @abstractmethod
    def scan_chunk(self, chunk: bytes, offset: int) -> list[TransitionHint]:
        """Scan a chunk for potential state transitions"""
        pass

    @abstractmethod
    def confirm_transition(
        self, data: bytes, hint: TransitionHint
    ) -> StateTransition | None:
        """Confirm a transition hint is valid and return exact positions"""
        pass
