from typing import Generic, TypeVar

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum

from atomea.io.text import StateScanner, StateTransition

S = TypeVar("S")


@dataclass
class ParsedRegion:
    """Container for parsed data from a state region"""

    state: Enum
    byte_range: tuple[int, int]
    data: dict


@dataclass
class ParsedFile:
    """Complete parsed file result"""

    file_path: str
    file_type: str
    regions: list[ParsedRegion]
    metadata: dict


class StateParser(ABC, Generic[S]):
    """Abstract parser for parsing specific state regions"""

    @abstractmethod
    def parse(self, data: bytes, state: S) -> dict:
        """Parse the bytes for a specific state region"""
        ...


class FileParser(ABC, Generic[S]):
    """Main parser orchestrator for a specific file type"""

    _scanner: StateScanner
    _parsers: dict[S, StateParser]

    def get_scanner(self) -> StateScanner:
        """Return the scanner for this file type"""
        return self._scanner

    def get_parser(self) -> StateParser[S] | None:
        """Return parser for specific state, None if no parsing needed"""
        return self._parsers

    def scan_bytes(self, buf: bytes) -> list[StateTransition]:
        transitions = self._scanner.scan(buf)
        return transitions

    def scan_file(self, file_path: str) -> list[StateTransition]:
        """Single pass through file to scan for states and their transitions."""
        with open(file_path, "rb") as fh:
            buf = fh.read()
        return self._scanner.scan(buf)

    def parse_file(self, file_path: str) -> ParsedFile:
        """Full parse of file"""
        ...


class ParserRegistry:
    """Registry for file type parsers"""

    def __init__(self):
        self._parsers: dict[str, FileParser] = {}

    def register(self, file_type: str, parser: FileParser) -> None:
        """Register a parser for a file type"""
        self._parsers[file_type] = parser

    def get_parser(self, file_type: str) -> FileParser | None:
        """Get parser for file type"""
        return self._parsers.get(file_type)
