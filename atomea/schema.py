from typing import Any

from collections.abc import Iterable

import yaml
from loguru import logger


class Atomea:
    """Define an atomistic schema."""

    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the schema.
        """
        self.schema: dict[str, Any] = {}
        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        if yaml_paths is not None:
            for yaml_path in yaml_paths:
                self.from_yaml(yaml_path)
        self.update(kwargs)
        self.validate()

    @property
    def schema(self) -> dict[str, Any]:
        return self._schema

    @schema.setter
    def schema(self, value: dict[str, Any]) -> None:
        if "description" in value.keys():
            value.pop("description")
        self._schema = value

    def from_yaml(self, yaml_path: str | None) -> None:
        """Load schema information from a YAML file. This will only update data
        contained in the YAML file.

        Args:
            yaml_path: Path to YAML file to load.
        """
        if yaml_path is not None:
            logger.info("Loading YAML schema from {}", yaml_path)
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.safe_load(f)
            logger.debug("YAML data:\n{}", yaml_data)
            self.update(yaml_data)
        self.yaml_path = yaml_path

    def update(self, attr_dict: dict[str, Any]) -> None:
        """Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        logger.debug("Updating schema:\n{}", attr_dict)
        for key, value in attr_dict.items():
            self.schema[key] = value

    def get(self) -> dict[str, Any]:
        """Retrieve the schema.

        Returns:
            A dictionary representing the current schema.
        """
        logger.debug("Retrieved schema:\n{}", self.schema)
        return self.schema

    def filter(self, keys: Iterable[str]) -> None:
        """Filter schema and only keep specified keys."""
        self.schema = {key: self.schema[key] for key in self.schema if key in keys}

    def fields(self, keys: Iterable[str] | None = None) -> dict[str, str]:
        """Return schema fields for tabular data."""
        schema = self.schema
        if keys is not None:
            schema = {key: schema[key] for key in schema if key in keys}
            return [(k, v["dtype"]) for k, v in schema.items()]
        return schema

    def validate(self) -> None:
        """Validate the schema."""
        schema = self.get()
        logger.debug("Validating schema:\n{}", schema)
        for key, value in schema.items():
            logger.trace("Validating {}", key)
            if isinstance(value, dict):
                for k in ("description", "ndim", "dtype", "units", "tabular", "length"):
                    if k not in value.keys():
                        raise ValueError(f"{key} is missing {k}")
            else:
                logger.warning(f"{key} is not a dictionary")

    def __enter__(self) -> dict[str, Any]:
        """Enter the context and return the current context as a dictionary."""
        return self.get()

    def __exit__(self, exc_type, exc_value, exc_tb):
        """Exit the context.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback information.
        """
