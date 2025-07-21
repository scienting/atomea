# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from enum import Enum, auto


class AmberV22State(Enum):
    """States for Amber v22 output files.

    This Enum defines the distinct logical sections or states that can be
    identified within an Amber 22 PMEMD output file. Each member represents
    a significant block of information, allowing for structured parsing.
    """

    UNKNOWN = auto()
    """An initial or undefined state."""

    PRELUDE = auto()
    """The introductory section of the file, typically containing
            version information (e.g., "Amber 22 PMEMD").
    """

    FILE_ASSIGNMENTS = auto()
    """Section detailing file input/output assignments."""

    INPUT_FILE = auto()
    """The section where the input file content is echoed."""

    SYSTEM_INFO = auto()
    """Information about resource usage and system details."""

    CONTROL_DATA = auto()
    """Parameters and control data used for the simulation run."""

    ATOMIC_COORDS = auto()
    """Section containing atomic coordinates and velocities."""

    RESULTS = auto()
    """
    The main results section, including energy values and other simulation outputs.
    """

    AVERAGES = auto()
    """Section presenting averaged values over the simulation."""

    TIMINGS = auto()
    """Information regarding the computational timings of the run."""

    COMPLETE = auto()
    """The final section, often indicating the completion of the run
    and total wall time.
    """


AMBER_V22_PATTERNS: dict[AmberV22State, list[bytes]] = {
    AmberV22State.PRELUDE: [
        rb"Amber 22 PMEMD",
    ],
    AmberV22State.FILE_ASSIGNMENTS: [rb"File Assignments:"],
    AmberV22State.INPUT_FILE: [rb"Here is the input file:"],
    AmberV22State.SYSTEM_INFO: [rb"   1.  RESOURCE   USE:"],
    AmberV22State.CONTROL_DATA: [rb"   2.  CONTROL  DATA  FOR  THE  RUN"],
    AmberV22State.ATOMIC_COORDS: [rb"   3.  ATOMIC COORDINATES AND VELOCITIES"],
    AmberV22State.RESULTS: [rb"   4.  RESULTS"],
    AmberV22State.AVERAGES: [rb"      A V E R A G E S   O V E R"],
    AmberV22State.TIMINGS: [rb"   5.  TIMINGS"],
    AmberV22State.COMPLETE: [rb"Master Total wall time:"],
}
"""A dictionary mapping `AmberV22State` members to lists of byte patterns.

These patterns are used by a `StateScanner` to detect transitions
between different states within an Amber 22 PMEMD output file.
Each pattern is a byte string (often a regex) that marks the beginning
of the corresponding state section.
"""
