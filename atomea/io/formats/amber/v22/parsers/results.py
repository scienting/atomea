# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

import numpy as np
from loguru import logger

from atomea.io.formats.amber.v22 import AmberV22State
from atomea.io.text import StateParser
from atomea.io.text.engines import ScanEngine, StdRegexEngine

PAT_NSTEP = rb"NSTEP\s*=\s*\d+"
PAT_ETOT = rb"^\s*Etot\s*="
PAT_BOND = rb"^\s*BOND\s*="
PAT_VOLUME = rb"\s*VOLUME\s*="
PAT_DENSITY = rb"^\s*Density\s*="


_PATTERN_MAP: dict[str, list[bytes]] = {
    "nstep": [PAT_NSTEP],
    "etot": [PAT_ETOT],
    "bond": [PAT_BOND],
    "volume": [PAT_VOLUME],
    "density": [PAT_DENSITY],
}


class AmberResultsParser(StateParser[AmberV22State]):
    """
    Parses an Amber `RESULTS` section using any ScanEngine implementation.
    """

    def __init__(self, engine: ScanEngine | None = None):
        """
        Args:
            engine: A concrete engine (e.g., StdRegexEngine).
                If omitted the std-lib regex engine is used.
        """
        self.engine: ScanEngine = engine or StdRegexEngine()
        self.engine.compile(_PATTERN_MAP)

    @classmethod
    def _split_num(cls, line: str) -> list[str]:
        return [tok for tok in line.split() if tok != "="]

    def parse(self, data: bytes, state: AmberV22State) -> dict[str, Any]:
        if state is not AmberV22State.RESULTS:
            logger.error("Invalid state for results parser")
            return {}

        # Temporary collectors (Python lists → transformed to NumPy at end)
        nsteps: list[int] = []
        times_ps: list[float] = []
        temps_k: list[float] = []
        press: list[float] = []
        etots: list[float] = []
        ektots: list[float] = []
        eptots: list[float] = []
        bonds: list[float] = []
        angles: list[float] = []
        diheds: list[float] = []
        volumes: list[float] = []
        densities: list[float] = []

        text = data.decode("utf-8", "ignore")  # decode once

        # Helper: quickly parse one *line* without another regex search

        matches = sorted(self.engine.stream(data), key=lambda t: t[0])

        # Build a dict pos→bucket so we can peek ahead
        pos_to_bucket = {start: bucket for start, _end, _, bucket in matches}
        line_starts = sorted(pos_to_bucket)

        for start in line_starts:
            bucket = pos_to_bucket[start]
            # slice bytes → line-text (cheap; region is not huge)
            line_end = data.find(b"\n", start)  # -1 ⇒ last line
            if line_end == -1:
                line_end = len(data)
            line = text[start:line_end]

            if bucket == "nstep":
                toks = self._split_num(line)
                # toks: ['NSTEP', '500', 'TIME(PS)', '1021.000', 'TEMP(K)', '302.32', 'PRESS', '0.0']
                nsteps.append(int(toks[1]))
                times_ps.append(float(toks[3]))
                temps_k.append(float(toks[5]))
                press.append(float(toks[7]))

            elif bucket == "etot":
                toks = self._split_num(line)
                etots.append(float(toks[1]))
                ektots.append(float(toks[3]))
                eptots.append(float(toks[5]))

            elif bucket == "bond":
                toks = self._split_num(line)
                bonds.append(float(toks[1]))
                angles.append(float(toks[3]))
                diheds.append(float(toks[5]))

            elif bucket == "volume":
                volumes.append(float(self._split_num(line)[1]))

            elif bucket == "density":
                densities.append(float(self._split_num(line)[1]))

        # Convert everything to arrays and return.
        if not nsteps:
            return {}

        result = {
            "nstep": np.array(nsteps),
            "time_ps": np.array(times_ps),
            "temp_k": np.array(temps_k),
            "pressure": np.array(press),
            "etot": np.array(etots),
            "ektot": np.array(ektots),
            "eptot": np.array(eptots),
            "bond": np.array(bonds),
            "angle": np.array(angles),
            "dihed": np.array(diheds),
            "volume": np.array(volumes),
            "density": np.array(densities),
            "n_steps": len(nsteps),
        }
        return result
