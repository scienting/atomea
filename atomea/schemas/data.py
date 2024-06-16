from typing import Any

from abc import ABC

import numpy as np


class Data(ABC):
    """
    Assists with updating data with atomistic schemas.
    """

    def update_atomistic(
        self,
        data: dict[str, Any],
        schema_map: dict[str, dict[str, str]],
        mol_index: int = 0,
    ) -> int:
        """
        Update the fields of the Schema instance with the provided data.

        This method updates the attributes of the Schema instance based on
        the keys and values in the provided dictionary. The keys in the dictionary
        can represent nested fields using dot notation.

        Args:
            data: A dictionary containing the keys and values to update the
                MoleculeSchema instance. The keys can use dot notation to
                specify nested attributes.
            schema_map: A mapping of field keys to their cadence and other metadata.
            mol_index: The current molecule index for updating array fields.

        Returns:
            The updated molecule index after processing the input data.

        Example:
            ```python
            mol_schema = MoleculeSchema()
            schema_map = mol_schema.get_schema_map()
            data = {
                "qc.energy": -76.4,
                "system.coordinates": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                "topology.bonds": [(0, 1), (0, 2)]
            }
            molecule.update(data, schema_map)
            ```

        Raises:
            AttributeError: If a specified attribute does not exist in the schema.
            ValueError: If the cadence is unknown.

        Notes:
            - The method supports updating nested attributes by splitting keys on the
              dot ('.') character.
            - If a key does not use dot notation, it will update the top-level
              attribute directly.
        """
        schema_map_alt = {v["field_key"]: v["cadence"] for k, v in schema_map.items()}
        for field_key, value in data.items():
            if field_key not in schema_map_alt.keys():
                continue

            cadence = schema_map_alt[field_key]
            if cadence == "molecule":
                mol_index = self._update_array(field_key, value, mol_index)
            elif cadence == "ensemble":
                self._set_field(field_key, value)
            else:
                raise ValueError(f"Unknown cadence: {cadence}")
        return mol_index

    def _set_field(self, key: str, value: Any, separator: str = ".") -> None:
        """
        Set a single field in the schema.

        Args:
            key: The key of the field to update. Can use dot notation for nested
                attributes.
            value: The value to set for the specified field.
            separator: The string used to separate nested keys.
        """
        keys = key.split(separator)
        if len(keys) > 1:
            sub_model = self
            for sub_key in keys[:-1]:
                sub_model = getattr(sub_model, sub_key)
            setattr(sub_model, keys[-1], value)
        else:
            setattr(self, key, value)

    def _get_field_value(self, key: str) -> Any:
        """
        Retrieve the value of a field in the schema.

        Args:
            key: The key of the field to retrieve. Can use dot notation for
                nested attributes.

        Returns:
            A tuple containing the keys, sub_model, and the current value of the field.
        """
        keys = key.split(".")
        sub_model = self
        for sub_key in keys[:-1]:
            sub_model = getattr(sub_model, sub_key)
        array = getattr(sub_model, keys[-1], None)
        return keys, sub_model, array

    def _update_array(self, key: str, value: Any, mol_index: int = 0) -> int:
        """
        Update an array field in the schema, resizing if necessary.

        Args:
            key: The key of the array field to update. Can use dot notation for
                nested attributes.
            value: The value to add to the array.
            mol_index: The current index for appending new values.

        Returns:
            The updated molecule index after appending the new values.

        Raises:
            TypeError: If the value is not a numpy array or cannot be converted to one.
        """
        keys, sub_model, array = self._get_field_value(key)

        if not isinstance(value, np.ndarray):
            value = np.array(value)
        if value.ndim == 2:
            value = value[None, ...]
        if array is None:
            array = np.empty((1000, *value.shape[1:]), dtype=value.dtype)

        # Check if we need to resize array to 2 times its current number of molecules.
        mol_index_new = int(mol_index + value.shape[0])
        if mol_index_new > array.shape[0]:
            new_shape = (int(array.shape[0] * 2), *array.shape[1:])
            new_array = np.empty(new_shape, dtype=array.dtype)
            new_array[:mol_index] = array
            array = new_array

        array[mol_index:mol_index_new] = value
        setattr(sub_model, keys[-1], array)

        return mol_index_new

    def _trim_molecule_arrays(
        self, mol_index: int, schema_map: dict[str, dict[str, str]]
    ) -> None:
        """
        Finalize arrays by trimming off the unused portions based on the current index.

        Args:
            mol_index: The current index up to which the arrays should be retained.
            schema_map: A mapping of field keys to their cadence and other metadata.
        """
        schema_map_alt = {v["field_key"]: v["cadence"] for k, v in schema_map.items()}
        for field_key, field_cadence in schema_map_alt.items():
            if field_cadence != "molecule":
                continue
            keys, sub_model, value = self._get_field_value(field_key)
            if not isinstance(value, np.ndarray):
                continue

            setattr(sub_model, keys[-1], value[:mol_index])
