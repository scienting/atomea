try:
    import MDAnalysis as mda

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
from typing import Any

from collections.abc import Collection

import numpy as np
import numpy.typing as npt

from ..schema import Atomea
from .digester import Digester


def accumate_things(*args, **kwargs):
    return list(*args)


class MDAnalysisDigester(Digester):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def digest(
        cls, atomea: Atomea, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        """Digest simulations supported by [MDAnalysis](https://www.mdanalysis.org/)."""
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is not installed")
        data: dict[str, Any] = {}
        u = mda.Universe(*args, **kwargs)
        with atomea as schema:
            for k in schema.keys():
                if k not in data.keys():
                    data[k] = []
                try:
                    data[k].append(getattr(cls, k)(u))
                except AttributeError as e:
                    raise NotImplementedError(f"{k} not implemented") from e
        return data

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
        atom_types: list[str] = u.atoms.accumulate("types", function=accumate_things)
        return atom_types
