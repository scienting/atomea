from abc import ABC


class Render(ABC):
    """Handles rendering files."""

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
