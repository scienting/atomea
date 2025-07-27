# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from atomea.io.formats.amber.v22 import (
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
    """A concrete `FileParser` implementation for Amber v22 output files.

    This parser specializes in understanding the structure of Amber 22 PMEMD
    output files. It uses a `StateScanner` configured with `AMBER_V22_PATTERNS`
    to identify different sections (states) within the file and then employs
    specific `StateParser` instances (e.g., `AmberSystemInfoParser`,
    `AmberResultsParser`) to extract structured data from those sections.

    Inherits from:
        FileParser[AmberV22State]: The generic base class for file parsers,
            specialized for `AmberV22State` enum.
    """

    def __init__(self):
        """Initializes the AmberV22Parser.

        Configures the internal `StateScanner` with Amber v22 specific patterns
        and sets up a dictionary of `StateParser` instances for handling
        different identified states within the Amber output file.
        """
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

        This method performs the complete parsing workflow for an Amber v22
        output file. It first reads the file, then scans it to identify
        state transitions, and finally iterates through these transitions
        to parse each relevant region using the appropriate `StateParser`.
        The extracted data is then encapsulated into `ParsedRegion` objects
        and aggregated into a `ParsedFile`.

        Args:
            file_path: Path to the Amber v22 output file to parse.

        Returns:
            A `ParsedFile` object containing all extracted data from the
                file, organized into `ParsedRegion` objects. This structured
                output is suitable for further processing or storage.
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
