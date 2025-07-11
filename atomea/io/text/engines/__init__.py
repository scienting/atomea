"""Engines for scanning and parsing bytes.

This package provides abstract interfaces and concrete implementations
for scanning byte streams to identify state transitions and parsing
regions of bytes into structured data.
"""

from .core import ScanEngine, S
from .std import StdRegexEngine

__all__ = ["ScanEngine", "S", "StdRegexEngine"]
