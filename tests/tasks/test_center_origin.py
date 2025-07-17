import os
import shutil

import numpy as np
import pytest
from raygent import TaskManager
from raygent.results.handlers import ResultsCollector

import atomea.typing as adt
from atomea.tasks.geom import CenterOriginTask, GeometricCenterTask


@pytest.fixture
def temp_project_dir(tmp_dir):
    """Provides a temporary directory for project stores."""
    dir_path = os.path.join(tmp_dir, "syn_geom_center")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    yield str(dir_path)


def test_generate_and_validate_valid_data(data_generator, temp_project_dir):
    """
    Test generating standard, valid synthetic data and validating it.
    """
    generator = data_generator(
        n_frames=(10, 50),
        n_atoms=(100, 100),
        random_seed=42,  # Ensure reproducibility for this test
    )
    batch_size = 10
    # Generate and write data to a new project
    project, reference_data = generator.create_synthetic_project(
        temp_project_dir,
        ens_id="test_ensemble",
        run_id="run_valid",
        clean_up_paths=True,
    )

    coords_ref: adt.Float64 = reference_data["ref_coords"]
    mean_ref = np.mean(coords_ref, axis=1)
    print(f"Reference coordinates: {coords_ref}")
    print(f"Reference mean coordinates: {mean_ref}")

    centers_ref = coords_ref - mean_ref[:, np.newaxis, :]
    print(f"Reference centers: {centers_ref}")

    manager = TaskManager[adt.Float64, ResultsCollector[adt.Float64]](
        CenterOriginTask, ResultsCollector, in_parallel=False
    )
    handler = manager.submit_tasks(
        project["test_ensemble"].coordinates.iter("run_valid", None, None, batch_size),
        batch_size=batch_size,
        prebatched=True,
    )
    results = handler.get()
    assert len(results) == 5
    assert np.allclose(results[0], centers_ref[:batch_size])
    assert np.allclose(results[1], centers_ref[batch_size : batch_size * 2])
