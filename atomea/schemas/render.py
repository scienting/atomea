from abc import ABC


class Render(ABC):
    """Handles rendering files."""

    dir_work: str = "."
    """
    Working directory to write files.
    """

    dir_input: str = "."
    """
    Path to a directory relative to
    [`dir_work`][schemas.RenderingConfig.dir_work] that will contain
    input files.
    """

    dir_output: str = "."
    """
    Path to a directory relative to
    [`dir_work`][schemas.RenderingConfig.dir_work] that the simulation will
    store output files.
    """

    def render(self) -> list[str]:
        """Prepare input lines by rendering templates or combining input configuration."""
        raise NotImplementedError

    def write_render(self, file_path: str) -> None:
        """Thin wrapper to write lines from the `render` function.

        Args:
            file_path: Path to write file.
        """
        lines = self.render()
        with open(file_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
