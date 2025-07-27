# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import pytest

from atomea.io.formats.amber.v22 import (
    AMBER_V22_PATTERNS,
    AmberV22Parser,
    AmberV22State,
    parsers,
)
from atomea.io.text import ParsedFile, StateScanner


class TestAmberV22Parser:
    def test_scan_complete_file(self, tmp_path, complete_amber_file):
        # Write test file
        test_file = tmp_path / "test_amber_scan.out"
        test_file.write_bytes(complete_amber_file)

        parser = AmberV22Parser()
        transitions = parser.scan_file(test_file)

        assert isinstance(transitions, list)

        # Check we found expected transitions
        states_found = [t.to_state for t in transitions]
        assert AmberV22State.PRELUDE in states_found
        assert AmberV22State.FILE_ASSIGNMENTS in states_found
        assert AmberV22State.RESULTS in states_found
        assert AmberV22State.TIMINGS in states_found

        # Check transitions are sorted by position
        positions = [t.start_pos for t in transitions]
        assert positions == sorted(positions)

        # Check end positions are filled in
        for i in range(len(transitions) - 1):
            assert transitions[i].end_pos == transitions[i + 1].start_pos

    def test_parse_complete_file(self, tmp_path, complete_amber_file):
        test_file = tmp_path / "test_amber_parse.out"
        test_file.write_bytes(complete_amber_file)

        parser = AmberV22Parser()
        parsed_file = parser.parse_file(test_file)

        assert isinstance(parsed_file, ParsedFile)

        assert parsed_file.file_type == "amber_v22"
        assert len(parsed_file.regions) == 2

        # Check metadata
        assert "n_regions" in parsed_file.metadata

        # Find energy regions
        results_regions = [
            r for r in parsed_file.regions if r.state == AmberV22State.RESULTS
        ]
        assert len(results_regions) == 1

        # Check first energy data
        first_energy = results_regions[0].data
        assert "nstep" in first_energy
        assert "etot" in first_energy

    def test_empty_file(self, tmp_path):
        empty_file = tmp_path / "empty.out"
        empty_file.write_bytes(b"")

        parser = AmberV22Parser()
        transitions = parser.scan_file(empty_file)

        assert isinstance(transitions, list)
        assert len(transitions) == 0

    def test_file_not_found(self):
        parser = AmberV22Parser()
        with pytest.raises(FileNotFoundError):
            parser.scan_file("/nonexistent/file.out")


class TestPatternMatching:
    """Test byte pattern matching performance and edge cases"""

    def test_pattern_at_chunk_boundary(self):
        """Test pattern split across chunk boundary"""
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)

        # Simulate pattern split across chunks
        chunk1 = b"some text   4.  RES"
        chunk2 = b"ULTS\n----------------"

        # Neither chunk alone should find the pattern
        hints1 = scanner.scan(chunk1)
        hints2 = scanner.scan(chunk2)

        # This is a limitation of simple scanning - would need overlap handling
        assert not any(h.pattern == b"   4.  RESULTS" for h in hints1)
        assert not any(h.pattern == b"   4.  RESULTS" for h in hints2)

    def test_case_sensitivity(self):
        """Ensure patterns are case-sensitive"""
        scanner = StateScanner(AMBER_V22_PATTERNS, AmberV22State.UNKNOWN)
        data = b"amber 22 pmemd"  # lowercase

        hints = scanner.scan(data)
        # Should not match because pattern is "Amber 22 PMEMD"
        assert not any(h.pattern == b"Amber 22 PMEMD" for h in hints)


@pytest.mark.parametrize(
    "nstep,time,temp,expected_nstep",
    [
        (0, 0.0, 300.0, 0),
        (999999, 1999.998, 298.15, 999999),
        (1, 0.002, 300.0, 1),
    ],
)
def test_energy_parser_parametrized(nstep, time, temp, expected_nstep):
    """Parametrized test for different energy values"""
    parser = parsers.AmberResultsParser()
    data = f" NSTEP = {nstep:8d}   TIME(PS) = {time:11.3f}  TEMP(K) = {temp:7.2f}  PRESS =     0.0\n".encode()

    result = parser.parse(data, AmberV22State.RESULTS)
    assert isinstance(result, dict)
    assert result["nstep"] == expected_nstep
    assert result["time_ps"] == time
    assert result["temp_k"] == temp
