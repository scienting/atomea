# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

"""Engines for scanning and parsing bytes.

This package provides abstract interfaces and concrete implementations
for scanning byte streams to identify state transitions and parsing
regions of bytes into structured data.
"""

from .core import ScanEngine, S
from .std import StdRegexEngine

__all__ = ["ScanEngine", "S", "StdRegexEngine"]
