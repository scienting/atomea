from dataclasses import dataclass
from enum import Enum

from atomea.digesters.text import StateParser


@dataclass
class ParseTask:
    """Task for parallel parsing"""

    state: Enum
    byte_range: tuple[int, int]
    parser: StateParser


@dataclass
class ParseResult:
    """Result from parallel parse task"""

    task: ParseTask
    result: dict
