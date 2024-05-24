from typing import Any

from collections.abc import Collection

import numpy as np
import numpy.typing as npt
from loguru import logger

try:
    import MDAnalysis as mda

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

from .digester import Digester


def accumulate_things(*args, **kwargs):
    """Helper function for accumulating data over multiple structures in a MDAnalysis
    Universe."""
    return list(*args)


class MDAnalysisDigester(Digester):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def checks():
        """Perform checks to ensure MDAnalysis is installed.

        Raises:
            ImportError: If MDAnalysis is not installed.
        """
        logger.info("Checking MDAnalysis installation")
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is not installed")

    @classmethod
    def prepare_step_inputs(
        cls, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, mda.Universe]:
        """Prepare and return the inputs necessary for the MDAnalysis digestion process.

        This method loads the MDAnalysis Universe using the provided arguments and keyword
        arguments and prepares the inputs for digestion.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            dict[str, mda.Universe]: A dictionary containing the MDAnalysis Universe.
        """
        logger.info("Preparing inputs for MDAnalysis digesters.")
        logger.debug(f"Loading universe with {args} and {kwargs}")
        if "coordinates" in kwargs.keys():
            logger.warning("Coordinates must be passed in args not kwargs")
            logger.warning(
                "Dropping coordinates from kwargs to avoid MDAnalysis error."
            )
            kwargs.pop("coordinates")
        u = mda.Universe(*args, **kwargs)
        logger.debug(f"MDAnalysis universe with {len(u.atoms)} atoms.")
        inputs = {"u": u}
        return inputs

    @staticmethod
    def coordinates(u: mda.Universe) -> dict[str, npt.NDArray[np.float64]]:
        """Return the coordinates of the atoms in the MDAnalysis Universe.

        Args:
            u (mda.Universe): The MDAnalysis Universe containing the atoms.

        Returns:
            dict[str, npt.NDArray[np.float64]]: A dictionary with the key
                "system.coordinates" and the value being a numpy array of atom
                coordinates.

        Raises:
            ValueError: If the coordinates are not a numpy array.
        """
        v = u.atoms.positions
        if isinstance(v, np.ndarray):
            return {"system.coordinates": v}
        else:
            raise ValueError("Coordinates must be a numpy array.")

    @staticmethod
    def ff_atom_type(u: mda.Universe) -> dict[str, list[str]]:
        """Return the force field atom types of the atoms in the MDAnalysis Universe.

        Args:
            u (mda.Universe): The MDAnalysis Universe containing the atoms.

        Returns:
            dict[str, list[str]]: A dictionary with the key "topology.ff_atom_type"
                and the value being a list of atom types.
        """
        atom_types: list[str] = u.atoms.accumulate("types", function=accumulate_things)
        return {"topology.ff_atom_type": atom_types}
