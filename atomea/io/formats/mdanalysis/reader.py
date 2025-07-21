# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import inspect

import numpy as np
from loguru import logger

try:
    import MDAnalysis as mda
    from MDAnalysis.guesser.tables import SYMB2Z

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

from typing import Any

from atomea.containers import Project
from atomea.helpers.bonding import get_molecule_ids
from atomea.io import Reader


class MDAnalysisReader(Reader):
    """
    Read data via MDAnalysis.
    """

    @classmethod
    def checks(cls) -> None:
        if not HAS_MDANALYSIS:
            raise ImportError("MDAnalysis is required for MDAnalysisReader")

    @classmethod
    def prepare(
        cls, topology: str, trajectory: str, *args: Any, **kwargs: Any
    ) -> dict[str, Any]:
        """Load a single Universe."""
        u = mda.Universe(topology, trajectory, *args, **kwargs)
        return {"u": u}

    @classmethod
    def extract(
        cls, prj: Project, ens_id: str, run_id: str, ctx: dict[str, Any]
    ) -> Project:
        """Extract and parse all possible information given the context.
        This method is responsible for calling all other methods.
        """
        for name, method in inspect.getmembers(cls, predicate=inspect.isfunction):
            if not name.startswith("parse_"):
                continue
            logger.debug(f"[{ens_id}] running parser `{name}`")
            try:
                prj = method(prj, ens_id, ctx, run_id=run_id)
            except Exception:
                logger.exception(f"[{ens_id}] error in parser `{name}`")
        return prj

    @staticmethod
    def parse_topology(
        prj: Project, ens_id: str, ctx: dict[str, Any], run_id: str | None = None
    ) -> Project:
        """Store topology."""
        u = ctx["u"]
        ids_component = np.array(u.atoms.resids, dtype=np.dtype(np.uint32))
        prj[ens_id].topology.ids.components.write(ids_component, run_id=run_id)

        labels_component = np.array(u.atoms.resnames, dtype=np.dtype(np.str_))
        prj[ens_id].topology.labels.components.write(labels_component, run_id=run_id)

        covalent = np.array(u.bonds.indices, dtype=np.dtype(np.uint64))
        prj[ens_id].topology.connectivity.bonds.covalent.write(covalent, run_id=run_id)

        angles = np.array(u.angles.indices, dtype=np.dtype(np.uint64))
        prj[ens_id].topology.connectivity.angles.write(angles, run_id=run_id)

        dihedrals = np.array(u.dihedrals.indices, dtype=np.dtype(np.uint64))
        prj[ens_id].topology.connectivity.dihedrals.write(dihedrals, run_id=run_id)

        if len(covalent) > 0:
            molecule_ids = get_molecule_ids(covalent, u.atoms.n_atoms)
            prj[ens_id].topology.ids.molecules.write(molecule_ids, run_id=run_id)

        return prj

    @staticmethod
    def parse_ensemble_metadata(
        prj: Project, ens_id: str, ctx: dict[str, Any], run_id: str | None = None
    ) -> Project:
        """Store per-ensemble atom symbols & ff types once."""
        u = ctx["u"]
        syms = np.array(u.atoms.elements, dtype=np.dtype(np.str_))
        prj[ens_id].topology.atoms.symbols.write(syms, run_id=run_id)
        types = np.array(u.atoms.types, dtype=np.dtype(np.str_))
        prj[ens_id].topology.atoms.types.write(types, run_id=run_id)
        return prj

    @staticmethod
    def parse_coordinates(
        prj: Project, ens_id: str, ctx: dict[str, Any], run_id: str | None = None
    ) -> Project:
        """Stack all frames into a (n_frames, n_atoms, 3) array."""
        u = ctx["u"]
        coords = [ts.positions.copy() for ts in u.trajectory]
        arr = np.stack(coords, axis=0, dtype=np.dtype(np.float64))
        # this writes a 3D microstate array
        prj[ens_id].coordinates.write(arr, run_id=run_id)
        return prj

    @staticmethod
    def parse_atom_z(
        prj: Project, ens_id: str, ctx: dict[str, Any], run_id: str | None = None
    ) -> Project:
        """Store the atomic numbers as an ensemble-level array."""
        u = ctx["u"]
        zs = np.array([SYMB2Z.get(sym, 0) for sym in u.atoms.elements])
        prj[ens_id].topology.atoms.atomic_numbers.write(zs, run_id=run_id)
        return prj
