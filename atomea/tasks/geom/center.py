from typing import Any, Type

import numpy as np
import numpy.typing as npt
from raygent.task import Task


class GeometricCenterTask(Task[npt.NDArray[np.float64], npt.NDArray[np.float64]]):
    """
    A Task to compute the geometric center for each structure (microstate)
    within a given chunk of atomic coordinates.

    This task is designed to work with chunks of coordinates, making it suitable
    for large datasets and parallel processing. It expects input as a NumPy array
    where each row represents a microstate and contains the flattened coordinates
    of all atoms for that microstate, or a 3D array (n_microstates, n_atoms, 3).

    The geometric center is calculated as the mean of the x, y, and z coordinates
    for all atoms in each structure.

    InputType: npt.NDArray[np.float64]
        A NumPy array representing a chunk of atomic coordinates.
        Expected shape: (n_microstates_in_chunk, n_atoms, 3)
        Alternatively, if flattened: (n_microstates_in_chunk, n_atoms * 3)

    OutputType: npt.NDArray[np.float64]
        A NumPy array where each row is the geometric center (x, y, z)
        for a corresponding microstate in the input chunk.
        Expected shape: (n_microstates_in_chunk, 3)

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
            at_once=True,  # Important: tells Task.run to call process_items
        )

        # Retrieve results
        geometric_centers = center_task_manager.get_results()
        ```
    """

    def process_items(
        self, items: npt.NDArray[np.float64], **kwargs: Any
    ) -> npt.NDArray[np.float64]:
        """
        Computes the geometric center for each structure in the input chunk.

        Args:
            items: Atomic coordinates of one or more microstates.
            **kwargs: Additional keyword arguments (not used in this specific task).

        Returns:
            A list of NumPy arrays, where each array is the geometric center (x, y, z)
            for a corresponding microstate. Each inner array will have shape (3,).
        """
        if items.ndim == 2:
            items = items[np.newaxis, :]
        if items.ndim != 3:
            raise RuntimeError("Items must be a NumPy array with 3 dimensions")
        centers = np.mean(items, axis=0)
        return centers
