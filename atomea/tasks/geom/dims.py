from typing import override

import numpy as np
from raygent import Task

import atomea.typing as adt


class BoxDimensionsTask(Task[adt.Float64, adt.Float64]):
    """
    A Task to compute the box dimensions of each structure (microstate)
    within a given chunk of atomic coordinates.


    """

    @override
    def do(self, batch: adt.Float64, *args: object, **kwargs: object) -> adt.Float64:
        """
        Computes the lengths of the structures assuming a box around
        the structure for each in the input chunk.

        Args:
            batch: Atomic coordinates of one or more microstates.
            **kwargs: Additional keyword arguments (not used in this specific task).
        Returns:
            A list of NumPy arrays, where each array is the box length (Lx, Ly, Lz)
            for a corresponding microstate. Each inner array will have shape (3,).
        """
        if batch.ndim == 2:
            batch = batch[np.newaxis]
        assert batch.ndim == 3, (
            "`batch` must be a NumPy array with 3 dimensions. Got {}",
            batch.ndim,
        )
        lengths = np.max(batch, axis=1) - np.min(batch, axis=1)
        return lengths
