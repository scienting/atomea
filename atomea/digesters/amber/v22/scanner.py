from atomea.digesters.amber.v22 import AMBER_V22_PATTERNS, AmberV22State
from atomea.digesters.text import (
    LiteralPattern,
    Pattern,
    StateScanner,
    StateTransition,
    TransitionHint,
)


class AmberV22Scanner(StateScanner[AmberV22State]):
    """Scanner for Amber v22 output files"""

    def __init__(self):
        self.patterns = AMBER_V22_PATTERNS
        self.current_state = AmberV22State.PRELUDE

    def get_patterns(self) -> dict[AmberV22State, list[Pattern]]:
        return self.patterns

    def scan_chunk(self, chunk: bytes, offset: int) -> list[TransitionHint]:
        """Scan chunk for potential state transitions"""
        hints = []

        # Simple implementation - scan for literal patterns
        for state, patterns in self.patterns.items():
            for pattern in patterns:
                if isinstance(pattern, LiteralPattern):
                    pos = chunk.find(pattern.bytes)
                    if pos != -1:
                        hints.append(
                            TransitionHint(
                                pattern_matched=pattern.bytes,
                                position=offset + pos,
                                confidence=1.0,
                                to_state=state,
                            )
                        )
        # Sort so we can get the transitions
        hints.sort(key=lambda t: t.position)
        for i in range(1, len(hints[1:])):
            hints[i].from_state = hints[i - 1].to_state

        return hints

    def confirm_transition(
        self, data: bytes, hint: TransitionHint
    ) -> StateTransition | None:
        """Confirm a transition hint and determine exact boundaries.

        Args:
            data: Content of chunk.
            hint:
        """
        # Find which state this pattern belongs to
        target_state = None
        for state, patterns in self.patterns.items():
            for pattern in patterns:
                if (
                    isinstance(pattern, LiteralPattern)
                    and pattern.bytes == hint.pattern_matched
                ):
                    target_state = state
                    break
            if target_state:
                break

        if not target_state:
            return None

        # Find the start of the line containing the pattern
        start_pos = hint.position
        while start_pos > 0 and data[start_pos - 1 : start_pos] != b"\n":
            start_pos -= 1

        # For section headers, we want to include the full header
        end_pos = hint.position + len(hint.pattern_matched)

        return StateTransition(
            from_state=self.current_state,
            to_state=target_state,
            start_pos=start_pos,
            end_pos=end_pos,
            trigger_pattern=hint.pattern_matched,
        )
