from atomea.io.amber.v22 import (
    AMBER_V22_PATTERNS,
    AmberV22State,
    parsers,
)
from atomea.io.text import (
    FileParser,
    ParsedFile,
    ParsedRegion,
    StateParser,
    StateScanner,
)


class AmberV22Parser(FileParser[AmberV22State]):
    def __init__(self):
        self._scanner = StateScanner(
            AMBER_V22_PATTERNS, first_state=AmberV22State.PRELUDE
        )
        """Scanner for state transitions"""
        self._parsers: dict[AmberV22State, StateParser] = {
            AmberV22State.SYSTEM_INFO: parsers.AmberSystemInfoParser(),
            AmberV22State.RESULTS: parsers.AmberResultsParser(),
        }
        """Parsers for regions identified by state."""

    def parse_file(self, file_path: str) -> ParsedFile:
        """Scan and parse file into regions with data ready to put
        into a store.

        Args:
            file_path: Path to file to parse.

        Returns:
            Parsed data from the file. You should use this to write
                into stores.
        """
        with open(file_path, "rb") as fh:
            buf = fh.read()
        transitions = self.scan_bytes(buf)

        regions = []
        for tr in transitions:
            p = self._parsers.get(tr.to_state)  # type: ignore
            if not p:
                continue
            data = buf[tr.start_pos : tr.end_pos]
            res = p.parse(data, tr.to_state)
            regions.append(
                ParsedRegion(
                    state=tr.to_state, byte_range=(tr.start_pos, tr.end_pos), data=res
                )
            )

        return ParsedFile(
            file_path=file_path,
            file_type="amber_v22",
            regions=regions,
            metadata={"n_regions": len(regions)},
        )
