from typing import Any

import inspect
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
from .ids import SchemaUUID


def accumulate_things(*args, **kwargs):
    """Helper function for accumulating data over multiple structures in a MDAnalysis
    Universe."""
    return list(*args)


class MDAnalysisDigester(Digester):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def get_uuid_map(cls) -> dict[str, str]:
        """
        Update the function UUID map by inspecting the class methods
        decorated with @SchemaUUID.
        """
        uuid_map = {}
        for name, method in inspect.getmembers(cls, predicate=inspect.isfunction):
            if name[:2] == "__":
                continue
            if callable(method) and hasattr(method, "__uuid__"):
                uuid_map[method.__uuid__] = name
        return uuid_map

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
    def prepare_inputs_digester(
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

    @classmethod
    def get_inputs_frame(cls, inputs_digester: dict[str, Any]) -> dict[str, Any]:
        inputs_frame = {"atoms": inputs_digester["u"].atoms}
        return inputs_frame

    @classmethod
    def next_frame(cls, inputs_digester: dict[str, Any]) -> dict[str, Any]:
        """Move the MDAnalysis Universe to the next trajectory frame.

        Args:
            inputs: A dictionary of inputs for the digestion process.

        Returns:
            A dictionary of inputs for the digestion process.
        """
        next(inputs_digester["u"].trajectory)
        return inputs_digester

    @staticmethod
    @SchemaUUID("81c7cec9-beec-4126-b6d8-91bee28951d6")
    def coordinates(atoms: mda.AtomGroup) -> npt.NDArray[np.float64]:
        """Return the coordinates of the atoms.

        Args:
            atoms: The MDAnalysis atoms object associated with the frame.

        Returns:
            A dictionary with the key `"system.coordinates"` and the value being a NumPy
            array of atom coordinates.

        Raises:
            ValueError: If the coordinates are not a numpy array.
        """
        v = atoms.positions
        if isinstance(v, np.ndarray):
            return v
        else:
            raise ValueError("Coordinates must be a numpy array.")

    @staticmethod
    @SchemaUUID("e34c0e1b-0eaa-4679-b060-3fcfe737aa15")
    def ff_atom_type(atoms: mda.AtomGroup) -> list[str]:
        """Return the force field atom types of the atoms.

        Args:
            atoms: The MDAnalysis atoms object associated with the frame.

        Returns:
            A dictionary with the key `"topology.ff_atom_type"` and the value being a
            list of atom types.
        """
        atom_types: list[str] = atoms.accumulate("types", function=accumulate_things)
        return atom_types
