from dataclasses import dataclass
from enum import Enum


@dataclass
class StateTransition:
    """Track byte positions of state changes for parsing."""

    from_state: Enum
    to_state: Enum
    start_pos: int
    end_pos: int
    pattern: bytes | None = None
