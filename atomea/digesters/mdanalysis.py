import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.guesser.tables import SYMB2Z

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

from typing import Any

from atomea.digesters.digester import Digester
from atomea.containers import Project


class MDAnalysisDigester(Digester):
    """
    Reads a trajectory via MDAnalysis; for each parse_… step,
    writes into the Project/Ensemble stores.
    """

    @classmethod
    def checks(cls) -> None:
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is required for MDAnalysisDigester")

    @classmethod
    def prepare_inputs(
        cls, topology: str, trajectory: str, **kwargs: Any
    ) -> dict[str, Any]:
        """Load a single Universe once."""
        u = mda.Universe(topology, trajectory, **kwargs)
        return {"u": u}

    @staticmethod
    def parse_ensemble_metadata(
        ctx: dict[str, Any], proj: Project, ens_id: str
    ) -> None:
        """Store per‐ensemble atom symbols & ff types once."""
        u = ctx["u"]
        syms = np.array(u.atoms.elements, dtype=str)
        types = np.array(u.atoms.types)
        proj.ensembles[ens_id].microstates.atom_symbol = syms
        proj.ensembles[ens_id].topology.ff_atom_type = types

    @staticmethod
    def parse_coordinates(ctx: dict[str, Any], proj: Project, ens_id: str) -> None:
        """Stack all frames into a (n_frames,n_atoms,3) array."""
        u = ctx["u"]
        coords = [ts.positions.copy() for ts in u.trajectory]
        arr = np.stack(coords, axis=0)
        # this writes a 3D microstate array
        proj.ensembles[ens_id].microstates.coordinates = arr

    @staticmethod
    def parse_atom_z(ctx: dict[str, Any], proj: Project, ens_id: str) -> None:
        """Store the atomic numbers as an ensemble‐level array."""
        u = ctx["u"]
        zs = np.array([SYMB2Z.get(sym, 0) for sym in u.atoms.elements])
        proj.ensembles[ens_id].microstates.atom_z = zs
