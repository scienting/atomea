from typing import Any

from loguru import logger
from pydantic import BaseModel

from ....io import YamlIO


class AmberSchemaBase(BaseModel, YamlIO):
    r"""Validate Amber contexts."""

    inputs: Any = NotImplemented

    def input_lines(self) -> list[str]:
        """
        Generate a list of input lines for the AMBER simulation.

        Returns:
            A list of strings representing the input lines.
        """
        input_lines = ["", "&cntrl"]

        for key, value in self.inputs.model_dump().items():
            if key in ("restraintmask", "timask1", "scmask1", "timask2", "scmask2"):
                value = f'"{value}"'
            line_to_add = f"    {key}={value},"
            logger.debug(f"Adding input line: {line_to_add.strip()}")
            input_lines.append(line_to_add)

        input_lines.append("&end")
        return input_lines

    def write_input(self, path: str) -> None:
        """
        Write the input file for the simulation.

        Args:
            path (str): The path to the output file.

        Returns:
            None
        """
        input_lines = self.input_lines()
        input_lines.append("")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(input_lines))
