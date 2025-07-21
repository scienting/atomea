# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

"""Processing AMBER simulations with roGFP2"""

import os
import shutil

import numpy as np

from atomea.containers import Project
from atomea.io.formats import MDAnalysisReader
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore


def test_reader_amber_rogfp2_serial(amber_rogfp2_sim_paths, tmp_dir):
    """
    Test the reading Amber simulations using the roGFP2 ensemble.

    Verifies that MDAnalysisReader builds a Project, creates a default Ensemble,
    and populates microstates and topology data correctly.
    """
    path_store_array = os.path.join(tmp_dir, "amber_rogfp2_serial.zarr")
    if os.path.exists(path_store_array):
        shutil.rmtree(path_store_array)

    path_store_table = os.path.join(tmp_dir, "amber_rogfp2_serial.tables")
    if os.path.exists(path_store_table):
        shutil.rmtree(path_store_table)

    store_array = ZarrArrayStore(path_store_array, mode="a")
    store_table = PolarsTableStore(
        path_store_table, mode="a", disk_format=DiskFormat.PARQUET
    )
    prj = Project(store_array, store_table)
    reader_args = (
        amber_rogfp2_sim_paths["mol.prmtop"],
        amber_rogfp2_sim_paths["07_relax_npt.nc"],
    )
    reader_kwargs = {
        "topology_format": "PRMTOP",
        "format": "NC",
    }
    prj: Project = MDAnalysisReader.run(
        prj,
        "default",
        "0",
        reader_args,  # type: ignore
        reader_kwargs,
    )
    ensemble = prj["default"]

    # coordinates: shape (n_frames, n_atoms, 3)
    coords = ensemble.coordinates.read(run_id="0")
    assert coords is not None
    assert coords.shape[0] == 100
    # spot‐check a few values
    assert np.allclose(coords[0, 0, 0], 33.924496)
    assert np.allclose(coords[0, 32, 0], 27.2496)
    assert np.allclose(coords[-1, 78, 0], 29.406982)

    # atom symbols
    syms = ensemble.topology.atoms.symbols.read(run_id="0")
    assert syms is not None
    assert syms[0] == "N"
    assert syms[8324] == "H"
    assert syms[-1] == "H"

    ids_component = ensemble.topology.ids.components.read(run_id="0")
    assert ids_component is not None
    assert ids_component[0] == 1
    assert ids_component[-1] == 10270

    labels_component = ensemble.topology.labels.components.read(run_id="0")
    assert labels_component is not None
    assert labels_component[0] == "MET"
    assert labels_component[-1] == "WAT"

    ids_molecules = ensemble.topology.ids.molecules.read(run_id="0")
    assert ids_molecules is not None
    assert ids_molecules[0] == 0
    assert ids_molecules[5729] == 701
    assert ids_molecules[-1] == 9985

    fftypes = ensemble.topology.atoms.types.read(run_id="0")
    assert fftypes is not None
    assert fftypes[0] == "N3"


def test_reader_write_amber_rogfp2_serial(amber_rogfp2_sim_paths, tmp_dir):
    """
    Test reader + storage into a Zarr store.
    """
    path_store_array = os.path.join(tmp_dir, "amber_rogfp2_serial.zarr")
    if os.path.exists(path_store_array):
        shutil.rmtree(path_store_array)

    path_store_table = os.path.join(tmp_dir, "amber_rogfp2_serial.tables")
    if os.path.exists(path_store_table):
        shutil.rmtree(path_store_table)

    store_array = ZarrArrayStore(path_store_array, mode="a")
    store_table = PolarsTableStore(
        path_store_table, mode="a", disk_format=DiskFormat.PARQUET
    )
    prj = Project(store_array, store_table)
    reader_args = (
        amber_rogfp2_sim_paths["mol.prmtop"],
        amber_rogfp2_sim_paths["07_relax_npt.nc"],
    )
    reader_kwargs = {
        "topology_format": "PRMTOP",
        "format": "NC",
    }
    MDAnalysisReader.run(
        prj,
        "default",
        "0",
        reader_args,  # type: ignore
        reader_kwargs,
    )

    # Re-open for reading and spot‐check
    store_r = ZarrArrayStore(path=path_store_array, mode="r")
    coords = store_r.read("default/0/coordinates")
    assert coords is not None
    assert np.allclose(coords[0, 0, 0], 33.924496)
    assert np.allclose(coords[0, 32, 0], 27.2496)
    assert np.allclose(coords[-1, 78, 0], 29.406982)
