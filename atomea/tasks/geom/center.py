from typing import override

import numpy as np
from raygent import Task

import atomea.typing as adt


class GeometricCenterTask(Task[adt.Float64, adt.Float64]):
    """
    A Task to compute the geometric center for each structure (microstate)
    within a given chunk of atomic coordinates.

    This task is designed to work with chunks of coordinates, making it suitable
    for large datasets and parallel processing. It expects input as a NumPy array
    where each row represents a microstate and contains the flattened coordinates
    of all atoms for that microstate, or a 3D array (n_microstates, n_atoms, 3).

    The geometric center is calculated as the mean of the x, y, and z coordinates
    for all atoms in each structure.

    InputType: adt.Float64
        A NumPy array representing a chunk of atomic coordinates.
        Expected shape: (n_microstates, n_atoms, 3)

    OutputType: adt.Float64
        A NumPy array where each row is the geometric center (x, y, z)
        for a corresponding microstate in the input chunk.
        Expected shape: (n_microstates, 3)

    Examples:

        ```python
        from raygent import TaskManager
        from atomea.tasks.geom import GeometricCenterTask

        # Instantiate the TaskManager with the new task
        center_task_manager = TaskManager(
            task_class=GeometricCenterTask,
            use_ray=False,
        )

        # Submit tasks using the .iter() method from Ensemble.coordinates
        # The chunk_size here determines how many microstates are processed at once
        # by the `process_items` method of GeometricCenterTask.
        center_task_manager.submit_tasks(
            items=ensemble.coordinates.iter(run_id="run_01", chunk_size=10),
            at_once=True,
        )

        # Retrieve results
        geometric_centers = center_task_manager.get_results()
        ```
    """

    @override
    def do(self, batch: adt.Float64, *args: object, **kwargs: object) -> adt.Float64:
        """
        Computes the geometric center for each structure in the input chunk.

        Args:
            items: Atomic coordinates of one or more microstates.
            **kwargs: Additional keyword arguments (not used in this specific task).

        Returns:
            A list of NumPy arrays, where each array is the geometric center (x, y, z)
            for a corresponding microstate. Each inner array will have shape (3,).
        """
        if batch.ndim == 2:
            batch = batch[np.newaxis, :]
        if batch.ndim != 3:
            raise RuntimeError("`batch` must be a NumPy array with 3 dimensions")
        centers = np.mean(batch, axis=1)
        return centers
