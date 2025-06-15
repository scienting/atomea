from enum import Enum, auto

from atomea.digesters.text import LiteralPattern, Pattern


class AmberV22State(Enum):
    """States for Amber v22 output files"""

    PRELUDE = auto()
    FILE_ASSIGNMENTS = auto()
    INPUT_FILE = auto()
    SYSTEM_INFO = auto()
    CONTROL_DATA = auto()
    ATOMIC_COORDS = auto()
    RESULTS = auto()
    AVERAGES = auto()
    TIMINGS = auto()
    COMPLETE = auto()


AMBER_V22_PATTERNS: dict[AmberV22State, list[Pattern]] = {
    AmberV22State.PRELUDE: [
        LiteralPattern(b"Amber 22 PMEMD"),
        LiteralPattern(b"-------------------------------------------------------"),
    ],
    AmberV22State.FILE_ASSIGNMENTS: [
        LiteralPattern(b"File Assignments:"),
    ],
    AmberV22State.INPUT_FILE: [
        LiteralPattern(b"Here is the input file:"),
    ],
    AmberV22State.SYSTEM_INFO: [
        LiteralPattern(b"   1.  RESOURCE   USE:"),
        LiteralPattern(b"getting box info from"),
    ],
    AmberV22State.CONTROL_DATA: [
        LiteralPattern(b"   2.  CONTROL  DATA  FOR  THE  RUN"),
    ],
    AmberV22State.ATOMIC_COORDS: [
        LiteralPattern(b"   3.  ATOMIC COORDINATES AND VELOCITIES"),
    ],
    AmberV22State.RESULTS: [
        LiteralPattern(b"   4.  RESULTS"),
    ],
    AmberV22State.AVERAGES: [
        LiteralPattern(b"      A V E R A G E S   O V E R"),
    ],
    AmberV22State.TIMINGS: [
        LiteralPattern(b"   5.  TIMINGS"),
    ],
    AmberV22State.COMPLETE: [
        LiteralPattern(b"Master Total wall time:"),
        LiteralPattern(b"Final Performance Info:"),
    ],
}
