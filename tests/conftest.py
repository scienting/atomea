import os
import shutil

import pytest
import numpy as np

from atomea import Project, enable_logging
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore

from .fixtures.amber_v22 import *  # type: ignore

TEST_DIR = os.path.dirname(__file__)
TMP_DIR = os.path.join(TEST_DIR, "tmp")
FILE_DIR = os.path.join(TEST_DIR, "files")


def pytest_sessionstart(session):
    r"""Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    # Creates a tmp directory for writing files.
    os.makedirs(TMP_DIR, exist_ok=True)


@pytest.fixture
def test_dir():
    return os.path.abspath(TEST_DIR)


@pytest.fixture
def tmp_dir():
    return os.path.abspath(TMP_DIR)


@pytest.fixture
def file_dir():
    return os.path.abspath(file_dir)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def amber_rogfp2_sim_paths():
    source_dir = "amber/rogfp2-sims"
    destination_dir = os.path.join(FILE_DIR, source_dir)

    paths = {
        "mol.prmtop": os.path.join(destination_dir, "mol.prmtop"),
        "07_relax_npt.nc": os.path.join(destination_dir, "07_relax_npt.nc"),
        "07_relax_npt.out": os.path.join(destination_dir, "07_relax_npt.out"),
    }
    return paths


@pytest.fixture
def test_project_with_synthetic_data():
    """
    Pytest fixture to create a Project with synthetic data for testing selectors.

    The synthetic system consists of:
    - Molecule 0: Water (OW, HW, HW) at origin, slightly shifted in microstate 1
    - Molecule 1: Water (OW, HW, HW) shifted +3 in x, slightly shifted in microstate 1
    - Molecule 2: Methane (C, H, H, H, H) shifted +6 in x, slightly shifted in microstate 1

    Total atoms: 11
    Total frames (microstates): 2
    """
    project_path = os.path.join(TMP_DIR, "test_project_selectors")
    os.makedirs(project_path, exist_ok=True)

    array_store_path = os.path.join(project_path, "syn.zarr")
    table_store_path = os.path.join(project_path, "syn.tables")

    # Ensure clean slate for stores
    if os.path.exists(array_store_path):
        shutil.rmtree(array_store_path, ignore_errors=True)
    if os.path.exists(table_store_path):
        shutil.rmtree(table_store_path, ignore_errors=True)

    array_store = ZarrArrayStore(array_store_path, mode="a")
    table_store = PolarsTableStore(
        table_store_path, mode="a", disk_format=DiskFormat.PARQUET
    )

    project = Project(array_store, table_store)
    ens_id = "test_ensemble"
    run_id = "test_run"
    ensemble = project.add_ensemble(ens_id)

    # Define synthetic data properties
    n_atoms = 11
    n_frames = 2

    # Atom Symbols, Types, Atomic Numbers (ensemble cadence data)
    atom_symbols = np.array(
        ["O", "H", "H", "O", "H", "H", "C", "H", "H", "H", "H"], dtype=np.dtype(np.str_)
    )
    atom_types = np.array(
        ["OW", "HW", "HW", "OW", "HW", "HW", "C", "H", "H", "H", "H"],
        dtype=np.dtype(np.str_),
    )
    atomic_numbers = np.array(
        [8, 1, 1, 8, 1, 1, 6, 1, 1, 1, 1], dtype=np.dtype(np.uint8)
    )

    # Molecule IDs, Component IDs, Component Labels (ensemble cadence data)
    mol_ids = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2], dtype=np.dtype(np.uint32))
    comp_ids = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2], dtype=np.dtype(np.uint32))
    comp_labels = np.array(
        [
            "Water",
            "Water",
            "Water",
            "Water",
            "Water",
            "Water",
            "Methane",
            "Methane",
            "Methane",
            "Methane",
            "Methane",
        ],
        dtype=np.dtype(np.str_),
    )

    # Coordinates (microstate cadence data)
    coords_ms0 = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.7, 0.7, 0.0],
            [0.7, -0.7, 0.0],  # Water 0 (O0, H1, H2)
            [3.0, 0.0, 0.0],
            [3.7, 0.7, 0.0],
            [3.7, -0.7, 0.0],  # Water 1 (O3, H4, H5)
            [6.0, 0.0, 0.0],
            [6.5, 0.5, 0.5],
            [6.5, -0.5, 0.5],
            [6.5, 0.5, -0.5],
            [6.5, -0.5, -0.5],  # Methane 2 (C6, H7, H8, H9, H10)
        ],
        dtype=np.dtype(np.float64),
    )

    coords_ms1 = coords_ms0 + np.array(
        [0.1, 0.0, 0.0]
    )  # All atoms shifted by +0.1 in x

    coordinates_data = np.stack([coords_ms0, coords_ms1])  # (n_frames, n_atoms, 3)

    # Write data to project
    ensemble.topology.atoms.symbols.write(run_id=run_id, data=atom_symbols)
    ensemble.topology.atoms.types.write(run_id=run_id, data=atom_types)
    ensemble.topology.atoms.atomic_numbers.write(run_id=run_id, data=atomic_numbers)
    ensemble.topology.ids.molecules.write(run_id=run_id, data=mol_ids)
    ensemble.topology.ids.components.write(run_id=run_id, data=comp_ids)
    ensemble.topology.labels.components.write(run_id=run_id, data=comp_labels)

    # Coordinates are written directly to the Data object as a stacked array
    ensemble.coordinates.write(run_id=run_id, data=coordinates_data)

    # Return the project, ensemble, run_id, n_atoms, and n_frames for tests to use
    yield project, ensemble, run_id, n_atoms, n_frames
