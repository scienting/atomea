import importlib
import os
import shutil

import pytest
import requests
from requests.compat import urljoin

from atomea import enable_logging
from atomea.schema import Atomea

TEST_DIR = os.path.dirname(__file__)
CACHE_DIR = os.path.join(TEST_DIR, "cache")


def download_file(url, local_path):
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    response = requests.get(url, timeout=10, stream=True)
    if response.status_code == 200:
        with open(local_path, "wb") as f:
            f.write(response.content)
    else:
        raise RuntimeError(f"Failed to download file at {url}")


def download_files(dest, base_url, slugs):
    paths = [os.path.join(dest, s) for s in slugs]
    if not os.path.exists(dest):
        os.makedirs(dest)
        urls = [urljoin(base_url, s) for s in slugs]
        try:
            for url, path in zip(urls, paths):
                download_file(url, path)
        except RuntimeError:
            shutil.rmtree(dest)
            raise
    return paths


@pytest.fixture
def test_dir():
    return os.path.abspath(TEST_DIR)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


@pytest.fixture
def base_atomea():
    schema_dir = importlib.resources.files("atomea.schemas")
    schema_paths = []
    for file_path in schema_dir.iterdir():
        if file_path.suffix == ".yaml" or file_path.suffix == ".yml":
            schema_paths.append(file_path)
    return Atomea(schema_paths)


@pytest.fixture
def globus_simlify_url():
    return "https://g-9e1ff7.1d26db.e229.dn.glob.us"


@pytest.fixture
def uuid_simlify_rogfp2(globus_simlify_url):  # pylint: disable=redefined-outer-name
    uuid = "f7498a8c-d021-491c-a343-10151e81434a"
    base_url = urljoin(globus_simlify_url, uuid) + "/"
    dest = os.path.join(CACHE_DIR, uuid)
    slugs = ["topo/mol.prmtop", "outputs/07_relax_npt.nc"]
    download_files(dest, base_url, slugs)
    return dest
