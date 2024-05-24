from typing import Any

import os
from collections.abc import Iterable, MutableSequence

import numpy as np

from ..schemas import EnsembleSchema
from ..utils import get_obj_from_string


class DiskData:
    """Manage data with a single backend."""

    def __init__(self, backend_array="zarr", backend_tabular="parquet"):
        self.backend_array = get_obj_from_string("atomea.io.backends._" + backend_array)
        self.backend_tabular = get_obj_from_string(
            "atomea.io.backends._" + backend_tabular
        )
        self.arrays = {}
        self.tables = {}
        self.array_keys: MutableSequence[str] = []
        self.table_atom_fields: dict[str, str] = {}
        self.table_structure_fields: dict[str, str] = {}

    def write_data(
        self, dest: str, idx: int, data: dict[str, Any], ensemble_schema: EnsembleSchema
    ) -> None:
        array_keys = ensemble_schema.storage_forms["array"]
        for k in array_keys:
            v = np.array(data[k]) if not isinstance(data[k], np.ndarray) else data[k]
            self.backend_array.write(self.arrays[k], v, idx)

        # We only write tables for the first structure. So we exit if this is not the
        # first structure.
        if idx > 0:
            return

        table_fields_atoms = ensemble_schema.table_fields["atoms"]
        if table_fields_atoms:
            data_path = os.path.join(dest, "atoms")
            data_tabular = {k: data[k] for k in table_fields_atoms.keys()}
            self.backend_tabular.write(data_path, data_tabular, table_fields_atoms)

        table_fields_structures = ensemble_schema.table_fields["structures"]
        if table_fields_structures:
            data_path = os.path.join(dest, "structures")
            data_tabular = {k: data[k] for k in table_fields_structures.keys()}
            self.backend_tabular.write(data_path, data_tabular, table_fields_structures)

    def store(
        self,
        dest: str,
        ensemble_schema: EnsembleSchema,
        digester: Any,
        digester_args: Iterable[Any],
        digester_kwargs: dict[str, Any],
    ) -> None:
        """Digest and store data as a stream with parallelization. Use this if you
        have a large amount of data to digest.

        Args:
            atomea: The atomea schema to use.
            digester: The digester to use.
            digester_args: Iterable of arguments to pass to the digester.
            digester_kwargs: Iterable of keyword arguments to pass to the digester.
        """
        schema = ensemble_schema.get()

        # Initialize all arrays
        self.array_keys = ensemble_schema.storage_forms["array"]
        self.arrays = {
            k: self.backend_array.initialize(
                dest,
                digester.array_size(schema[k], *digester_args, **digester_kwargs),
                schema[k]["dtype"],
            )
            for k in self.array_keys
        }

        # Note, we do not need to analyze tables because there will only be one
        # for "atoms" and "structures".

        for idx, data in enumerate(
            digester.digest_step(ensemble_schema, *digester_args, **digester_kwargs)
        ):
            self.write_data(dest, idx, data, ensemble_schema)
