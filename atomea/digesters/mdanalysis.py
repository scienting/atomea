from typing import Any

from collections.abc import Collection

import numpy as np
import numpy.typing as npt

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
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is not installed")

    @classmethod
    def prepare_step_inputs(
        cls, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, mda.Universe]:
        inputs = {"u": mda.Universe(*args, **kwargs)}
        return inputs

    @staticmethod
    def coordinates(u: mda.Universe) -> npt.NDArray[np.float64]:
        """Return the coordinates of the atoms in the universe."""
        v = u.atoms.positions
        if isinstance(v, np.ndarray):
            return v
        raise TypeError(f"{type(v)} is not a numpy array")

    @staticmethod
    def ff_atom_type(u: mda.Universe) -> list[str]:
        """Return the coordinates of the atoms in the universe."""
        atom_types: list[str] = u.atoms.accumulate("types", function=accumulate_things)
        return atom_types
