from dataclasses import dataclass
from enum import Enum


@dataclass
class StateTransition:
    """Track exact byte positions of state changes"""

    from_state: Enum
    to_state: Enum
    start_pos: int
    end_pos: int
    trigger_pattern: bytes | None = None


@dataclass
class TransitionHint:
    """Potential state transition found during scanning"""

    pattern_matched: bytes
    position: int
    confidence: float
