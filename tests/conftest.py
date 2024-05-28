import os

import pytest
from google.cloud import storage

from atomea import enable_logging

TEST_DIR = os.path.dirname(__file__)
TMP_DIR = os.path.join(TEST_DIR, "tmp")
CACHE_DIR = os.path.join(TEST_DIR, "cache")
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = os.getenv(
    "GOOGLE_APPLICATION_CREDENTIALS", ""
)


def pytest_sessionstart(session):  # pytest_configure(config)
    r"""Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    # Creates a tmp directory for writing files.
    os.makedirs(TMP_DIR, exist_ok=True)


def download_from_gcs(bucket_name, source_blob_name, destination_file_name):
    """Download a file from GCS bucket."""
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)


@pytest.fixture
def test_dir():
    return os.path.abspath(TEST_DIR)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def amber_rogfp2_sim_paths():
    bucket_name = "atomea-tests"
    source_dir = "amber_rogfp2_sims"
    destination_dir = os.path.join(CACHE_DIR, source_dir)

    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    paths = {
        "mol.prmtop": os.path.join(destination_dir, "mol.prmtop"),
        "07_relax_npt.nc": os.path.join(destination_dir, "07_relax_npt.nc"),
    }
    for file_name, save_path in paths.items():
        source_blob_name = os.path.join(source_dir, file_name)
        if not os.path.exists(save_path):
            download_from_gcs(bucket_name, source_blob_name, save_path)
    return paths
