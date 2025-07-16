import os
import shutil

import numpy as np

import pytest
from raygent import TaskManager
from raygent.results.handlers import ResultsCollector

import atomea.typing as adt
from atomea.tasks.geom import GeometricCenterTask


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
        n_frames=(5, 10),
        n_atoms=(20, 50),
        random_seed=42,  # Ensure reproducibility for this test
    )

    # Generate and write data to a new project
    project, reference_data = generator.create_synthetic_project(
        temp_project_dir,
        ens_id="test_ensemble",
        run_id="run_valid",
        clean_up_paths=True,
    )

    coords_ref: adt.Float64 = reference_data["ref_coords"]
    mean_ref = np.mean(coords_ref, axis=1)

    manager = TaskManager[adt.Float64, ResultsCollector[adt.Float64]](
        GeometricCenterTask, ResultsCollector, in_parallel=False
    )
    handler = manager.submit_tasks(
        project["test_ensemble"].coordinates.iter("run_valid", None, None, 5),
        batch_size=5,
        prebatched=True,
    )
    results = handler.get()

    assert len(results) == 2
    assert np.allclose(results[0], mean_ref[:5])
    assert np.allclose(results[1], mean_ref[5:])
    
