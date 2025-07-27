# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import numpy as np
import pytest

from atomea.io.formats.amber.v22 import AmberV22State, parsers


class TestAmberResultsParser:
    def test_parse_results(self, results_step_data):
        parser = parsers.AmberResultsParser()
        result = parser.parse(results_step_data, AmberV22State.RESULTS)

        assert len(result) > 0

        # Check parsed values
        assert np.isclose(result["nstep"][0], 500)
        assert np.isclose(result["nstep"][-1], 2000)
        assert np.isclose(result["time_ps"][0], 1021.0)
        assert np.isclose(result["time_ps"][-1], 1024.000)
        assert np.isclose(result["temp_k"][0], 302.32)
        assert np.isclose(result["temp_k"][-1], 300.74)
        assert np.isclose(result["pressure"][0], 0.0)
        assert np.isclose(result["pressure"][-1], 0.0)
        assert np.isclose(result["etot"][0], -109458.5000)
        assert np.isclose(result["etot"][-1], -109711.6013)
        assert np.isclose(result["ektot"][0], 20735.8122)
        assert np.isclose(result["ektot"][-1], 20627.5611)
        assert np.isclose(result["eptot"][0], -130194.3121)
        assert np.isclose(result["eptot"][-1], -130339.1625)
        assert np.isclose(result["bond"][0], 678.7182)
        assert np.isclose(result["bond"][-1], 669.1943)
        assert np.isclose(result["angle"][0], 1926.4968)
        assert np.isclose(result["angle"][-1], 1912.8495)
        assert np.isclose(result["dihed"][0], 1263.9976)
        assert np.isclose(result["dihed"][-1], 1239.1669)
        assert np.isclose(result["volume"][0], 332147.5960)
        assert np.isclose(result["volume"][-1], 332147.5960)
        assert np.isclose(result["density"][0], 1.0361)
        assert np.isclose(result["density"][-1], 1.0361)

    def test_parse_wrong_state(self, results_step_data):
        parser = parsers.AmberResultsParser()
        result = parser.parse(results_step_data, AmberV22State.SYSTEM_INFO)
        assert len(result) == 0

    def test_parse_incomplete_data(self):
        parser = parsers.AmberResultsParser()
        incomplete_data = b" NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   302.32  PRESS =     0.0n"

        with pytest.raises(ValueError):
            parser.parse(incomplete_data, AmberV22State.RESULTS)

    def test_parse_malformed_data(self):
        parser = parsers.AmberResultsParser()
        malformed_data = b"This is not energy data at all"

        result = parser.parse(malformed_data, AmberV22State.RESULTS)
        assert len(result) == 0


class TestAmberSystemInfoParser:
    def test_parse_system_info(self, system_info_data):
        parser = parsers.AmberSystemInfoParser()
        result = parser.parse(system_info_data, AmberV22State.SYSTEM_INFO)

        assert isinstance(result, dict)

        assert result["natom"] == 33582
        assert result["nres"] == 10270
        # Box dimensions not in this snippet
        assert "box" not in result
