from typing import Callable

from dataclasses import dataclass


@dataclass
class Pattern:
    """Base pattern type"""

    pass


@dataclass
class LiteralPattern(Pattern):
    bytes: bytes


@dataclass
class RegexPattern(Pattern):
    pattern: str


@dataclass
class CustomPattern(Pattern):
    matcher: Callable[[bytes], bool]
