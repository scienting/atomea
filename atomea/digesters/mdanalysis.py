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
from .uuids import SchemaUUID


def accumulate_things(*args, **kwargs):
    """Helper function for accumulating data over multiple structures in a MDAnalysis
    Universe."""
    return list(*args)


class MDAnalysisDigester(Digester):
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
            A dictionary containing the MDAnalysis Universe.
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
        """Return the atomic coordinates in Angstroms.

        **Schema UUID:** [`81c7cec9-beec-4126-b6d8-91bee28951d6`]
        [schemas.atomistic.system.SystemSchema.coordinates]

        Args:
            atoms: The MDAnalysis atoms object associated with the frame.

        Returns:
            A 2D array containing the Cartesian coordinates of all atoms in.

        Raises:
            ValueError: If the coordinates are not a 2D NumPy array.
        """
        v = atoms.positions
        if isinstance(v, np.ndarray):
            if v.ndim != 2:
                raise ValueError("Coordinates are not two dimensions")
            return v
        else:
            raise TypeError("Coordinates must be a numpy array.")

    @staticmethod
    @SchemaUUID("81c21a83-4b72-48c6-a576-4541b468eb90")
    def atom_symbols(atoms: mda.AtomGroup) -> npt.NDArray[np.str_]:
        """Return the atomic coordinates in Angstroms.

        **Schema UUID:** [`81c21a83-4b72-48c6-a576-4541b468eb90`]
        [schemas.atomistic.system.SystemSchema.atom_symbol]

        Args:
            atoms: The MDAnalysis atoms object associated with the frame.

        Returns:
            A 1D array containing the atomic symbols.
        """
        v = atoms.elements
        if isinstance(v, np.ndarray):
            return v.astype(np.str_)
        else:
            raise TypeError("atom_symbols must be a numpy array.")

    @staticmethod
    @SchemaUUID("e34c0e1b-0eaa-4679-b060-3fcfe737aa15")
    def ff_atom_type(atoms: mda.AtomGroup) -> list[str]:
        """Return the force field atom types of the atoms.

        **Schema UUID:** [`e34c0e1b-0eaa-4679-b060-3fcfe737aa15`]
        [schemas.atomistic.topology.TopologySchema.ff_atom_type]

        Args:
            atoms: The MDAnalysis atoms object associated with the frame.

        Returns:
            Force field atom types.
        """
        atom_types: list[str] = atoms.accumulate("types", function=accumulate_things)
        return atom_types
