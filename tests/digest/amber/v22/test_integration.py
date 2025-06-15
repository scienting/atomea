import pytest

from atomea.digesters.amber.v22 import (
    AmberV22Parser,
    AmberV22Scanner,
    AmberV22State,
    parsers,
)
from atomea.digesters.text import ParsedFile


class TestAmberV22Parser:
    def test_scan_complete_file(self, tmp_path, complete_amber_file):
        # Write test file
        test_file = tmp_path / "test_amber.out"
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
        test_file = tmp_path / "test_amber.out"
        test_file.write_bytes(complete_amber_file)

        parser = AmberV22Parser()
        parsed_file = parser.parse_file(test_file)

        assert isinstance(parsed_file, ParsedFile)

        assert parsed_file.file_type == "amber_v22"
        assert len(parsed_file.regions) == 2

        # Check metadata
        assert "file_type" in parsed_file.metadata
        assert "n_transitions" in parsed_file.metadata
        assert "states_found" in parsed_file.metadata

        # Find energy regions
        energy_regions = [
            r for r in parsed_file.regions if r.state == AmberV22State.RESULTS
        ]
        assert len(energy_regions) >= 2  # Should find both MD steps

        # Check first energy data
        if energy_regions:
            first_energy = energy_regions[0].data
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
        with pytest.raises(ValueError):
            parser.scan_file("/nonexistent/file.out")


class TestPatternMatching:
    """Test byte pattern matching performance and edge cases"""

    def test_pattern_at_chunk_boundary(self):
        """Test pattern split across chunk boundary"""
        scanner = AmberV22Scanner()

        # Simulate pattern split across chunks
        chunk1 = b"some text   4.  RES"
        chunk2 = b"ULTS\n----------------"

        # Neither chunk alone should find the pattern
        hints1 = scanner.scan_chunk(chunk1, 0)
        hints2 = scanner.scan_chunk(chunk2, len(chunk1))

        # This is a limitation of simple scanning - would need overlap handling
        assert not any(h.pattern_matched == b"   4.  RESULTS" for h in hints1)
        assert not any(h.pattern_matched == b"   4.  RESULTS" for h in hints2)

    def test_multiple_identical_patterns(self):
        """Test handling of repeated patterns"""
        scanner = AmberV22Scanner()
        data = b""" NSTEP =      500
 NSTEP =     1000
 NSTEP =     1500"""

        hints = scanner.scan_chunk(data, 0)
        nstep_hints = [h for h in hints if h.pattern_matched == b" NSTEP ="]

        # Should find only the first occurrence with simple find()
        assert len(nstep_hints) == 1

    def test_case_sensitivity(self):
        """Ensure patterns are case-sensitive"""
        scanner = AmberV22Scanner()
        data = b"amber 22 pmemd"  # lowercase

        hints = scanner.scan_chunk(data, 0)
        # Should not match because pattern is "Amber 22 PMEMD"
        assert not any(h.pattern_matched == b"Amber 22 PMEMD" for h in hints)


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
    assert result.value["nstep"] == expected_nstep
    assert result.value["time_ps"] == time
    assert result.value["temp_k"] == temp
