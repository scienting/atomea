import os

from atomea.schemas.simulation.amber import Amber22Schema

from .conftest import TMP_DIR


def test_render_amber_write():
    schema_amber = Amber22Schema()
    yaml_path = os.path.join(TMP_DIR, "amber22.yml")
    schema_amber.to_yaml(yaml_path)
