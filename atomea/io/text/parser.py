from typing import Any, Generic, TypeVar

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum

from atomea.io.text import StateScanner, StateTransition

S = TypeVar("S")


@dataclass
class ParsedRegion:
    """Container for parsed data from a state region"""

    state: Enum
    """
    The Enum member representing the state of the content from which this region
    was parsed.
    """

    byte_range: tuple[int, int]
    """
    A tuple `(start_byte, end_byte)` indicating the inclusive start and exclusive
    end byte positions of this region within the original file.
    """

    data: dict[str, Any]
    """
    A dictionary containing the parsed data specific to this
    region's state. The structure and content of this dictionary
    will vary depending on the `StateParser` used for this region.
    """


@dataclass
class ParsedFile:
    """Complete parsed file result.

    This dataclass encapsulates all information extracted from a file
    after a full parsing operation, including its path, type, a list
    of parsed regions, and any global metadata.
    """

    file_path: str
    """The absolute or relative path to the file that was parsed."""

    file_type: str
    """A string identifier for the type of file that was parsed (e.g., "amber_v22")."""

    regions: list[ParsedRegion]
    """
    A list of `ParsedRegion` objects, each representing a
    distinct, parsed section of the file. The order of regions
    in this list corresponds to their order of appearance in the file.
    """

    metadata: dict[str, Any]
    """
    A dictionary containing any high-level metadata about
    the parsed file, such as the number of regions found.
    """


class StateParser(ABC, Generic[S]):
    """Abstract parser for parsing specific state regions.

    This abstract base class defines the interface for parsers that are
    responsible for extracting structured data from a specific byte
    range corresponding to a particular state within a file.

    Type variables:
        S: An Enum type representing the possible states a file can be in.
    """

    @abstractmethod
    def parse(self, data: bytes, state: S) -> dict[str, Any]:
        """Parse the bytes for a specific state region.

        This method should be implemented by concrete `StateParser`
        subclasses to transform raw bytes from a file region into a
        structured dictionary representation.

        Args:
            data: The bytes corresponding to the specific region of the file
                that needs to be parsed.
            state: The `Enum` member indicating the current state of the
                content, which guides how the `data` should be parsed.

        Returns:
            A dictionary containing the structured data extracted from the
                provided `data` bytes, relevant to the given `state`.
        """


class FileParser(ABC, Generic[S]):
    """Main parser orchestrator for a specific file type.

    This abstract base class provides the core functionality for scanning
    a file to identify state transitions and then parsing the identified
    regions using appropriate `StateParser` instances.

    Type variables:
        S: An Enum type representing the possible states a file can be in.
    """

    _scanner: StateScanner
    """An instance of `StateScanner` used to identify state transitions within the file.
    """

    _parsers: dict[S, StateParser]
    """A dictionary mapping `Enum` states to their corresponding
    `StateParser` instances, responsible for parsing data within those states.
    """

    def get_scanner(self) -> StateScanner:
        """Return the scanner for this file type.

        Returns:
            The `StateScanner` instance configured for this file type.
        """
        return self._scanner

    def get_parser(self) -> StateParser[S] | None:
        """Return parser for specific state, None if no parsing needed.

        This method is intended to return the collection of state parsers,
        not a single parser for a specific state. It might be a placeholder
        or an incomplete implementation.

        Returns:
            The dictionary of `StateParser` instances, or `None` if no
                parsers are configured (though the current implementation
                always returns the dictionary).
        """
        return self._parsers

    def scan_bytes(self, buf: bytes) -> list[StateTransition]:
        """Scan a byte buffer to identify state transitions.

        This method delegates the scanning operation to the internal
        `StateScanner` to find all state changes within a given byte
        sequence.

        Args:
            buf: The byte buffer to be scanned.

        Returns:
            A list of `StateTransition` objects, detailing the start and
                end positions of each state and the transition patterns.
        """
        transitions = self._scanner.scan(buf)
        return transitions

    def scan_file(self, file_path: str) -> list[StateTransition]:
        """Single pass through file to scan for states and their transitions.

        Reads the entire file into memory and then scans it for state
        transitions.

        Args:
            file_path: The path to the file to be scanned.

        Returns:
            A list of `StateTransition` objects, detailing the start and
            end positions of each state and the transition patterns.
        """
        with open(file_path, "rb") as fh:
            buf = fh.read()
        return self._scanner.scan(buf)

    @abstractmethod
    def parse_file(self, file_path: str) -> ParsedFile:
        """Full parse of file.

        This abstract method must be implemented by concrete `FileParser`
        subclasses to orchestrate the complete parsing process: scanning
        the file, identifying regions, and then parsing each region using
        the appropriate `StateParser`.

        Args:
            file_path: The path to the file to be fully parsed.

        Returns:
            A `ParsedFile` object containing all extracted data and metadata.
        """
