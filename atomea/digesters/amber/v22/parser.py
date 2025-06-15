from loguru import logger

from atomea.digesters.amber.v22 import AmberV22Scanner, AmberV22State, parsers
from atomea.digesters.text import (
    FileParser,
    ParsedFile,
    ParsedRegion,
    StateParser,
    StateScanner,
    StateTransition,
)


class AmberV22Parser(FileParser[AmberV22State]):
    """Main parser for Amber v22 output files"""

    def __init__(self):
        self.scanner = AmberV22Scanner()
        self.parsers = {
            AmberV22State.MD_STEP: parsers.AmberResultsParser(),
            AmberV22State.SYSTEM_INFO: parsers.AmberSystemInfoParser(),
        }

    def get_scanner(self) -> StateScanner[AmberV22State]:
        return self.scanner

    def get_parser(self, state: AmberV22State) -> StateParser[AmberV22State] | None:
        return self.parsers.get(state)

    def scan_file(self, file_path: str) -> list[StateTransition]:
        """Scan file and return all state transitions"""
        try:
            transitions = []

            with open(file_path, "rb") as f:
                file_data = f.read()

            # Scan entire file
            hints = self.scanner.scan_chunk(file_data, 0)

            # Confirm each hint
            for hint in hints:
                transition = self.scanner.confirm_transition(file_data, hint)
                if transition:
                    transitions.append(transition)

            # Sort by position
            transitions.sort(key=lambda t: t.start_pos)

            # Fill in end positions
            for i in range(len(transitions) - 1):
                transitions[i].end_pos = transitions[i + 1].start_pos

            if transitions:
                transitions[-1].end_pos = len(file_data)

            return transitions

        except Exception as e:
            raise ValueError(f"Failed to scan file: {str(e)}")

    def parse_file(self, file_path: str) -> ParsedFile:
        """Full parse of Amber output file"""
        # First scan for transitions
        transitions = self.scan_file(file_path)

        regions = []

        try:
            with open(file_path, "rb") as f:
                file_data = f.read()

            # Parse each region
            for transition in transitions:
                parser = self.get_parser(transition.to_state)
                if parser:
                    # Extract region data
                    region_data = file_data[transition.start_pos : transition.end_pos]
                    parse_result = parser.parse(region_data, transition.to_state)

                    if parse_result:
                        regions.append(
                            ParsedRegion(
                                state=transition.to_state,
                                byte_range=(transition.start_pos, transition.end_pos),
                                data=parse_result,
                            )
                        )
                    else:
                        logger.warning("Failed to parse region")

            # Extract metadata
            metadata = {
                "file_type": "amber_v22",
                "n_transitions": len(transitions),
                "states_found": list(set(t.to_state for t in transitions)),
            }

            return ParsedFile(
                file_path=file_path,
                file_type="amber_v22",
                regions=regions,
                metadata=metadata,
            )

        except Exception as e:
            raise ValueError(f"Failed to parse file: {str(e)}")
