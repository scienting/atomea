from enum import Enum, auto


class AmberV22State(Enum):
    """States for Amber v22 output files"""

    UNKNOWN = auto()
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
