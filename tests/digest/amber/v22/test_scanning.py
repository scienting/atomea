from atomea.digesters.amber.v22 import AmberV22Scanner, AmberV22State
from atomea.digesters.text import TransitionHint


class TestAmberV22Scanner:
    def test_scanner_initialization(self):
        scanner = AmberV22Scanner()
        assert scanner.current_state == AmberV22State.PRELUDE
        assert len(scanner.patterns) > 0

    def test_find_amber_header(self, amber_header):
        scanner = AmberV22Scanner()
        hints = scanner.scan_chunk(amber_header, 0)

        # Should find the Amber 22 PMEMD pattern
        assert len(hints) >= 1
        assert any(hint.pattern_matched == b"Amber 22 PMEMD" for hint in hints)

    def test_find_results_step(self, results_step_data):
        scanner = AmberV22Scanner()
        hints = scanner.scan_chunk(results_step_data, 0)

        # Should find NSTEP pattern
        assert len(hints) >= 1
        assert any(hint.pattern_matched == b"   4.  RESULTS" for hint in hints)

    def test_find_multiple_patterns(self, complete_amber_file):
        scanner = AmberV22Scanner()
        hints = scanner.scan_chunk(complete_amber_file, 0)

        # Should find multiple patterns
        pattern_texts = [hint.pattern_matched for hint in hints]
        assert b"Amber 22 PMEMD" in pattern_texts
        assert b"File Assignments:" in pattern_texts
        assert b"   4.  RESULTS" in pattern_texts
        assert b"   5.  TIMINGS" in pattern_texts

    def test_confirm_transition(self, results_step_data):
        scanner = AmberV22Scanner()
        hint = TransitionHint(
            pattern_matched=b"   4.  RESULTS",
            position=results_step_data.find(b"   4.  RESULTS"),
            confidence=1.0,
        )

        transition = scanner.confirm_transition(results_step_data, hint)
        assert transition is not None
        assert transition.to_state == AmberV22State.RESULTS
        assert transition.trigger_pattern == b"   4.  RESULTS"
        # Should capture from start of line
        assert (
            results_step_data[transition.start_pos : transition.start_pos + 14]
            == b"   4.  RESULTS"
        )

    def test_scan_chunk_with_offset(self):
        scanner = AmberV22Scanner()
        chunk = b"Some prefix text\n    4.  RESULTS"
        offset = 1000  # Simulating chunk from middle of file

        hints = scanner.scan_chunk(chunk, offset)
        assert len(hints) == 1
        assert hints[0].position == offset + chunk.find(b"   4.  RESULTS")
