from typing import Generic

from atomea.digesters.text import StateTransition
from atomea.digesters.text.engines import S, ScanEngine, StdRegexEngine


class StateScanner(Generic[S]):
    """Single-pass scanner for state transitions."""

    def __init__(
        self,
        patterns: dict[S, list[bytes]],
        first_state: S,
        engine: ScanEngine[S] | None = None,
    ):
        self.first_state = first_state
        """Initial state of the file."""
        if engine is None:
            engine = StdRegexEngine()
        self.engine = engine
        self.engine.compile(patterns)

    @property
    def patterns(self):
        return self.engine._patterns

    def scan(self, buf: bytes) -> list[StateTransition]:
        """Single-pass scan of bytes to identify states and their transitions."""
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
