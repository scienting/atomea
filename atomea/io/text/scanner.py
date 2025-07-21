# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Generic

from atomea.io.text import StateTransition
from atomea.io.text.engines import S, ScanEngine, StdRegexEngine


class StateScanner(Generic[S]):
    """Single-pass scanner for state transitions.

    This class orchestrates the scanning of a byte buffer or file
    to identify transitions between different states based on predefined
    patterns. It uses a `ScanEngine` to perform the actual pattern matching.

    Type variables:
        S: An Enum type representing the possible states a file can be in.
    """

    def __init__(
        self,
        patterns: dict[S, list[bytes]],
        first_state: S,
        engine: ScanEngine[S] | None = None,
    ):
        """Initializes the StateScanner.

        Args:
            patterns: A dictionary where keys are state Enum members and
                values are lists of byte patterns (regex or literal) that
                signify a transition *to* that state.
            first_state: The initial state of the file, which is used to
                begin the scanning process.
            engine: An optional `ScanEngine` instance to use for pattern
                matching. If `None`, `StdRegexEngine` will be used by default.
        """
        self.first_state = first_state
        if engine is None:
            engine = StdRegexEngine()
        self.engine = engine
        self.engine.compile(patterns)

    @property
    def patterns(self):
        """Returns the compiled patterns used by the internal scan engine.

        This property provides access to the patterns that the scanner
        is configured to detect.

        Returns:
            The underlying patterns from the scan engine.
        """
        return self.engine._patterns

    def scan(self, buf: bytes) -> list[StateTransition]:
        """Single-pass scan of bytes to identify states and their transitions.

        This method iterates through the provided byte buffer, identifying
        all occurrences of predefined patterns that indicate a change in state.
        It constructs a list of `StateTransition` objects, each detailing
        the `from_state`, `to_state`, and byte `start_pos`/`end_pos` of
        each identified region.

        Args:
            buf: The byte buffer to be scanned for state transitions.

        Returns:
            A list of `StateTransition` objects, ordered by their appearance
                in the buffer, representing the sequence of states and their
                corresponding byte ranges.
        """
        cur_state, trans, last = self.first_state, [], 0

        for start, _end, pattern, next_state in sorted(
            self.engine.stream(buf), key=lambda t: t[0]
        ):
            trans.append(StateTransition(cur_state, next_state, start, -1, pattern))
            cur_state, last = next_state, start

        for i in range(len(trans) - 1):
            trans[i].end_pos = trans[i + 1].start_pos
        if trans:
            trans[-1].end_pos = len(buf)

        return trans
