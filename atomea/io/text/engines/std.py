# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import re
from collections.abc import Iterator

from atomea.io.text.engines import S, ScanEngine


class StdRegexEngine(ScanEngine[S]):
    """Engine that uses the Python standard library's regex engine.

    This `ScanEngine` implementation leverages Python's `re` module
    to perform pattern matching. It compiles regular expressions from
    byte patterns and uses them to find state transitions within text
    decoded from byte buffers.

    Type variables:
        S: An Enum type representing the possible states a file can be in.
    """

    def __init__(self):
        """Initializes the StdRegexEngine.

        Sets up internal lists to store uncompiled and compiled regex patterns
        along with their associated states.
        """
        self._patterns: list[tuple[re.Pattern, S]] = []
        self._compiled: list[tuple[re.Pattern, S]] = []

    def compile(self, patterns: dict[S, list[bytes]]) -> None:
        """Compile regex patterns to speed up searching.

        This method takes a dictionary of state-to-pattern mappings, decodes
        the byte patterns into strings, and compiles them into `re.Pattern`
        objects. These compiled patterns are stored internally for efficient
        streaming.

        Args:
            patterns: A dictionary mapping an enum representing the states
                a file could be in to a list of byte patterns (regex) that
                signify a transition to that state.
        """
        self._patterns = patterns
        for state, pat_list in patterns.items():
            for pat in pat_list:
                self._compiled.append((re.compile(pat.decode(), re.MULTILINE), state))

    def stream(self, buf: bytes) -> Iterator[tuple[int, int, bytes, S]]:
        """Stream (match_start, match_end, pattern, target_state) triples.

        Decodes the input byte buffer into a UTF-8 string (ignoring errors)
        and then uses the compiled regular expressions to find all matches.
        For each match, it yields the start and end positions, the original
        byte pattern that matched, and the associated target state.

        Args:
            buf: The byte buffer to be scanned.

        Yields:
            A tuple containing:

                - `int`: The inclusive starting byte offset of the match.
                - `int`: The exclusive ending byte offset of the match.
                - `bytes`: The original byte pattern that was matched.
                - `S`: The target state associated with the matched pattern.
        """
        text = buf.decode("utf-8", "ignore")
        for rex, state in self._compiled:
            for m in rex.finditer(text):
                yield m.start(), m.end(), m.re.pattern.encode(), state
