from typing import Any

import re

import numpy as np
from loguru import logger

from atomea.digesters.amber.v22 import AmberV22State
from atomea.digesters.text import StateParser


class AmberResultsParser(StateParser[AmberV22State]):
    """Parser for MD results steps"""

    def parse(self, data: bytes, state: AmberV22State) -> dict[str, Any]:
        """Parse all MD steps from the RESULTS section"""
        result: dict[str, Any] = {}

        if state != AmberV22State.RESULTS:
            logger.error("Invalid state for results parser")
            return result

        try:
            text = data.decode("utf-8")

            # Find all MD steps in the results section
            nstep_pattern = re.compile(
                r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([\d.]+)\s+TEMP\(K\)\s*=\s*([\d.]+)\s+PRESS\s*=\s*([\d.-]+)"
            )
            energy_pattern = re.compile(
                r"Etot\s*=\s*([\d.-]+)\s+EKtot\s*=\s*([\d.-]+)\s+EPtot\s*=\s*([\d.-]+)"
            )
            bond_pattern = re.compile(
                r"BOND\s*=\s*([\d.-]+)\s+ANGLE\s*=\s*([\d.-]+)\s+DIHED\s*=\s*([\d.-]+)"
            )
            volume_pattern = re.compile(r"VOLUME\s*=\s*([\d.-]+)")
            density_pattern = re.compile(r"Density\s*=\s*([\d.-]+)")

            # Split into individual step blocks (each starts with NSTEP)
            step_blocks = re.split(r"(?=\sNSTEP\s*=)", text)
            step_blocks = [block for block in step_blocks if block.strip()]

            # Initialize lists to collect data
            nsteps = []
            times = []
            temps = []
            pressures = []
            etots = []
            ektots = []
            eptots = []
            bonds = []
            angles = []
            diheds = []
            volumes = []
            densities = []

            # Parse each step block
            for block in step_blocks:
                # Parse NSTEP line
                nstep_match = nstep_pattern.search(block)
                if nstep_match:
                    nsteps.append(int(nstep_match.group(1)))
                    times.append(float(nstep_match.group(2)))
                    temps.append(float(nstep_match.group(3)))
                    pressures.append(float(nstep_match.group(4)))
                else:
                    continue  # Skip blocks without NSTEP

                # Parse energy components
                energy_match = energy_pattern.search(block)
                if energy_match:
                    etots.append(float(energy_match.group(1)))
                    ektots.append(float(energy_match.group(2)))
                    eptots.append(float(energy_match.group(3)))
                else:
                    # Use NaN for missing values
                    etots.append(np.nan)
                    ektots.append(np.nan)
                    eptots.append(np.nan)

                # Parse bond, angle, dihedral
                bond_match = bond_pattern.search(block)
                if bond_match:
                    bonds.append(float(bond_match.group(1)))
                    angles.append(float(bond_match.group(2)))
                    diheds.append(float(bond_match.group(3)))
                else:
                    bonds.append(np.nan)
                    angles.append(np.nan)
                    diheds.append(np.nan)

                # Parse volume and density
                volume_match = volume_pattern.search(block)
                if volume_match:
                    volumes.append(float(volume_match.group(1)))
                else:
                    volumes.append(np.nan)

                density_match = density_pattern.search(block)
                if density_match:
                    densities.append(float(density_match.group(1)))
                else:
                    densities.append(np.nan)

            # Convert lists to numpy arrays
            if nsteps:
                result["nstep"] = np.array(nsteps, dtype=np.int32)
                result["time_ps"] = np.array(times, dtype=np.float64)
                result["temp_k"] = np.array(temps, dtype=np.float64)
                result["pressure"] = np.array(pressures, dtype=np.float64)
                result["etot"] = np.array(etots, dtype=np.float64)
                result["ektot"] = np.array(ektots, dtype=np.float64)
                result["eptot"] = np.array(eptots, dtype=np.float64)
                result["bond"] = np.array(bonds, dtype=np.float64)
                result["angle"] = np.array(angles, dtype=np.float64)
                result["dihed"] = np.array(diheds, dtype=np.float64)
                result["volume"] = np.array(volumes, dtype=np.float64)
                result["density"] = np.array(densities, dtype=np.float64)

                # Add metadata
                result["n_steps"] = len(nsteps)

            return result

        except Exception as e:
            raise RuntimeError(f"Failed to parse energy data: {str(e)}")
