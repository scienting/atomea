# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from dataclasses import dataclass
from enum import Enum


@dataclass
class StateTransition:
    """Track byte positions of state changes for parsing.

    This dataclass represents a transition from one state to another within
    a byte stream or file, marking the boundaries of a specific content region.
    """

    from_state: Enum
    """
    The `Enum` member representing the state of the content immediately preceding
    this transition.
    """

    to_state: Enum
    """
    The `Enum` member representing the state the content is transitioning to,
    which applies to the region starting at `start_pos`.
    """

    start_pos: int
    """The inclusive byte offset in the original data where this new state
    region begins.
    """

    end_pos: int
    """The exclusive byte offset in the original data where this state region ends.
    This marks the beginning of the next state region or the end of the file.
    """

    pattern: bytes | None = None
    """
    An optional byte string representing the pattern that triggered this
    state transition. This can be useful for debugging or understanding the
    parsing logic.
    """
