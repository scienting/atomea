from typing import Any

from loguru import logger

from atomea.containers import Project
from atomea.digesters import Digester
from atomea.digesters.text import (
    FileParser,
    ParsedFile,
    ParsedRegion,
    StateParser,
    StateScanner,
)


class AmberOutputDigester(Digester):
    @classmethod
    def checks(cls) -> None:
        pass

    @classmethod
    def prepare(
        cls, parser: FileParser, *args: Any, **kwargs: dict[str, Any]
    ) -> dict[str, Any]:
        return {"parser": parser}

    @classmethod
    def extract(cls, prj: Project, id_ens: str, ctx: dict[str, Any]) -> Project:
        """Extract and parse all possible information given the context.
        This method is responsible for calling all other methods.
        """
        parser = ctx["parser"]
        parsed_file = parser.parse_file
        try:
            prj = method(prj, id_ens, ctx)
        except Exception:
            logger.exception(f"[{id_ens}] error in parser `{name}`")
        return prj

    @staticmethod
    def store_results_sections(prj: Project, region: ParsedRegion) -> Project:
        prj.energy.potential_mm.append(region.data["eptot"])
        return prj
