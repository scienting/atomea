# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from atomea.io.formats.amber.v22 import AMBER_V22_PATTERNS, AmberV22State
from atomea.io.text import StateScanner


class TestAmberV22Scanner:
    def test_scanner_initialization(self):
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)
        assert scanner.first_state == AmberV22State.UNKNOWN
        assert len(scanner.patterns) > 0

    def test_find_amber_header(self, amber_header):
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)
        transitions = scanner.scan(amber_header)

        # Should find the Amber 22 PMEMD pattern
        assert len(transitions) >= 1
        assert any(trans.pattern == b"Amber 22 PMEMD" for trans in transitions)

    def test_find_results_step(self, results_step_data):
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)
        transitions = scanner.scan(results_step_data)

        # Should find NSTEP pattern
        assert len(transitions) == 1
        assert any(trans.pattern == b"   4.  RESULTS" for trans in transitions)

    def test_find_multiple_patterns(self, complete_amber_file):
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)
        transitions = scanner.scan(complete_amber_file)

        # Should find multiple patterns
        pattern_texts = [trans.pattern for trans in transitions]
        assert b"Amber 22 PMEMD" in pattern_texts
        assert b"File Assignments:" in pattern_texts
        assert b"   4.  RESULTS" in pattern_texts
        assert b"   5.  TIMINGS" in pattern_texts
