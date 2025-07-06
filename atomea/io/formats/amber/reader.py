from typing import Any

from atomea.containers import Project
from atomea.io import Reader
from atomea.io.text import (
    FileParser,
    ParsedRegion,
)


class AmberOutputReader(Reader):
    @classmethod
    def checks(cls) -> None:
        pass

    @classmethod
    def prepare(
        cls, file_path: str, parser: FileParser, *args: Any, **kwargs: Any
    ) -> dict[str, Any]:
        return {"file_path": file_path, "parser": parser}

    @classmethod
    def extract(
        cls, prj: Project, ens_id: str, run_id: str, ctx: dict[str, Any]
    ) -> Project:
        """Extract and parse all possible information given the context.
        This method is responsible for calling all other methods.

        Args:
            prt: Project to store extracted data to.
            ens_id: Ensemble to append extracted data to.
            run_id: Run ID to append extracted data to.
            ctx: Context information for parsing data from.
        """
        parser = ctx["parser"]
        parsed_file = parser.parse_file(ctx["file_path"])

        microstate_id_next = -1  # placeholder until we determine

        for region in parsed_file.regions:
            for key, value in region.data.items():
                if key == "eptot":
                    data_obj = prj.energy.potential_mm
                    microstate_id_next = data_obj.next_microstate_id(ens_id, run_id)
                    df_next = data_obj.prep_dataframe(
                        ens_id, run_id, microstate_id_next, value
                    )
                    if microstate_id_next == 0:
                        data_obj.write(df_next, run_id=run_id)
                    else:
                        data_obj.append(df_next, run_id=run_id)

        return prj

    @staticmethod
    def store_results_sections(
        prj: Project, region: ParsedRegion, run_id: str | None = None
    ) -> Project:
        prj.energy.potential_mm.append(region.data["eptot"], run_id=run_id)
        return prj
