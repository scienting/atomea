from typing import Generic, TypeVar

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum

from atomea.digesters.text import StateScanner, StateTransition

T = TypeVar("T")


class StateParser(ABC, Generic[T]):
    """Abstract parser for parsing specific state regions"""

    @abstractmethod
    def parse(self, data: bytes, state: T) -> dict:
        """Parse the bytes for a specific state region"""
        pass


class FileParser(ABC, Generic[T]):
    """Main parser orchestrator for a specific file type"""

    @abstractmethod
    def get_scanner(self) -> StateScanner[T]:
        """Return the scanner for this file type"""
        pass

    @abstractmethod
    def get_parser(self, state: T) -> StateParser[T] | None:
        """Return parser for specific state, None if no parsing needed"""
        pass

    def scan_file(self, file_path: str) -> list[StateTransition]:
        """Scan entire file and return all state transitions"""
        ...

    def parse_file(self, file_path: str) -> dict:
        """Full parse of file"""
        ...


# Parser registry for managing multiple file types
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
