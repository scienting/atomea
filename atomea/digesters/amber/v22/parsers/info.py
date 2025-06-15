from typing import Any

import re

from atomea.digesters.amber.v22 import AmberV22State
from atomea.digesters.text import StateParser


class AmberSystemInfoParser(StateParser[AmberV22State]):
    """Parser for system information section"""

    def parse(self, data: bytes, state: AmberV22State) -> dict[str, Any]:
        """Parse system information"""
        if state != AmberV22State.SYSTEM_INFO:
            raise ValueError("Invalid state for system info parser")

        try:
            text = data.decode("utf-8")
            result: dict[str, Any] = {}

            # Parse atom counts
            natom_match = re.search(r"NATOM\s*=\s*(\d+)", text)
            if natom_match:
                result["natom"] = int(natom_match.group(1))

            # Parse residue count
            nres_match = re.search(r"NRES\s*=\s*(\d+)", text)
            if nres_match:
                result["nres"] = int(nres_match.group(1))

            # Parse box dimensions
            box_match = re.search(
                r"Box X\s*=\s*([\d.]+)\s+Box Y\s*=\s*([\d.]+)\s+Box Z\s*=\s*([\d.]+)",
                text,
            )
            if box_match:
                result["box"] = {
                    "x": float(box_match.group(1)),
                    "y": float(box_match.group(2)),
                    "z": float(box_match.group(3)),
                }

            # Parse water count
            water_match = re.search(
                r"Number of triangulated 3-point waters found:\s*(\d+)", text
            )
            if water_match:
                result["n_waters"] = int(water_match.group(1))

            return result

        except Exception as e:
            raise RuntimeError(f"Failed to parse system info: {str(e)}") from e
