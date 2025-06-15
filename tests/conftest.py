import os

import pytest

from atomea import enable_logging

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
    }
    return paths
