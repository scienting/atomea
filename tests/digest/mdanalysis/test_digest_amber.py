"""Processing AMBER simulations with roGFP2"""

import os
import shutil

import numpy as np

from atomea.containers import Project
from atomea.io import MDAnalysisReader
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
    store_table = PolarsTableStore(path_store_table, DiskFormat.PARQUET)
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
    ensemble = prj.ensembles["default"]

    # coordinates: shape (n_frames, n_atoms, 3)
    coords = ensemble.microstates.coordinates.value
    assert coords is not None
    assert coords.shape[0] == 100
    # spot‐check a few values
    assert np.allclose(coords[0, 0, 0], 33.924496)
    assert np.allclose(coords[0, 32, 0], 27.2496)
    assert np.allclose(coords[-1, 78, 0], 29.406982)

    # atom symbols
    syms = ensemble.microstates.atom_symbol.value
    assert syms is not None
    assert syms[0] == "N"
    assert syms[8324] == "H"
    assert syms[-1] == "H"

    # force-field atom types in topology
    fftypes = ensemble.topology.ff_atom_type.value
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
    store_table = PolarsTableStore(path_store_table, DiskFormat.PARQUET)
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
    coords = store_r.read("default/microstate/coordinates")
    assert coords is not None
    assert np.allclose(coords[0, 0, 0], 33.924496)
    assert np.allclose(coords[0, 32, 0], 27.2496)
    assert np.allclose(coords[-1, 78, 0], 29.406982)
