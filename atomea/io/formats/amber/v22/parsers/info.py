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

PAT_NATOM = rb"NATOM\s*=\s*\d+"
PAT_NRES = rb"NRES\s*=\s*\d+"
PAT_WATER = rb"Number of triangulated 3-point waters found:\s*\d+"

_PATTERN_MAP: dict[str, list[bytes]] = {
    "natom": [PAT_NATOM],
    "nres": [PAT_NRES],
    "water": [PAT_WATER],
}


def _split_num(line: str) -> list[str]:
    return [tok for tok in line.split() if tok != "="]


# ── 2.  Parser class ────────────────────────────────────────────────────────
class AmberSystemInfoParser(StateParser[AmberV22State]):
    """
    Parses the ‘SYSTEM INFO’ section (Amber v22) with a pluggable ScanEngine.
    """

    def __init__(self, engine: ScanEngine | None = None) -> None:
        # default = pure-Python std-lib regex backend
        self.engine: ScanEngine = engine or StdRegexEngine()
        # compile once per instance
        self.engine.compile(_PATTERN_MAP)

    # ----------------------------------------------------------------------
    def parse(self, data: bytes, state: AmberV22State) -> dict[str, Any]:
        if state is not AmberV22State.SYSTEM_INFO:
            logger.error("Invalid state for system-info parser")
            return {}

        text = data.decode("utf-8", "ignore")

        result: dict[str, Any] = {
            "natom": np.nan,
            "nres": np.nan,
            "n_water": np.nan,
        }

        # stream matches once, already bucket-labelled
        for start, _end, _match, bucket in self.engine.stream(data):
            # grab the whole line
            line_end = data.find(b"\n", start)
            if line_end == -1:
                line_end = len(data)
            line = text[start:line_end]

            if bucket == "natom":
                result["natom"] = int(_split_num(line)[1])

            elif bucket == "nres":
                result["nres"] = int(_split_num(line)[1])

            elif bucket == "water":
                result["n_water"] = int(_split_num(line)[-1])

        # remove items that remained NaN / not found
        clean = {
            k: v
            for k, v in result.items()
            if not (isinstance(v, float) and np.isnan(v))
        }

        return clean
