from typing import Any

import yaml


class YamlIO:
    def update(self, data: dict[str, Any]) -> None:
        """Update the fields of the Schema instance with the provided data.

        This method updates the attributes of the Schema instance based on
        the keys and values in the provided dictionary. The keys in the dictionary
        can represent nested fields using dot notation.

        Args:
            data (dict[str, Any]): A dictionary containing the keys and values to
            update the MoleculeSchema instance. The keys can use dot notation to
            specify nested attributes.

        Example:
            >>> molecule = MoleculeSchema()
            >>> update_data = {
            >>>     "identification.name": "Water",
            >>>     "qc.energy": -76.4,
            >>>     "system.coordinates": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            >>>     "topology.bonds": [(0, 1), (0, 2)]
            >>> }
            >>> molecule.update(update_data)

        Raises:
            AttributeError: If a specified attribute does not exist in the schema.

        Notes:
            -   The method supports updating nested attributes by splitting keys on the
                dot ('.') character.
            -   If a key does not use dot notation, it will update the top-level
                attribute directly.
        """
        for key, value in data.items():
            keys = key.split(".")
            if len(keys) > 1:
                sub_model = self
                for sub_key in keys[:-1]:
                    sub_model = getattr(sub_model, sub_key)
                setattr(sub_model, keys[-1], value)
            else:
                setattr(self, key, value)

    def from_yaml(self, yaml_paths: str | list[str]) -> None:
        """Update the instance's attributes from one or more YAML files.

        Args:
            yaml_paths: A sequence of YAML file paths or a single YAML file path.

        Raises:
            FileNotFoundError: If any of the YAML files cannot be found.
        """
        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        for yaml_path in yaml_paths:
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.safe_load(f)
                for key, value in yaml_data.items():
                    if key in self.__fields__:  # type: ignore # pylint: disable=no-member
                        setattr(self, key, value)

    def to_yaml(self, file_path: str) -> None:
        """Serialize a Pydantic BaseModel instance to a YAML file.

        Args:
            file_path: Path to the YAML file to write the serialized data to.

        Raises:
            IOError: If the file cannot be written to.
        """
        config_dict = self.dict()  # type: ignore # pylint: disable=no-member
        with open(file_path, "w", encoding="utf-8") as f:
            yaml.dump(config_dict, f, default_flow_style=False)
