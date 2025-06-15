from typing import Iterator

import re

from atomea.digesters.text.engines import S, ScanEngine


class StdRegexEngine(ScanEngine[S]):
    """Engine that uses the Python standard library's regex engine."""

    def __init__(self):
        self._patterns: list[tuple[re.Pattern, S]] = []
        self._compiled: list[tuple[re.Pattern, S]] = []

    def compile(self, patterns: dict[S, list[bytes]]) -> None:
        """Compile regex patterns to speed up searching.

        This sets the `_compiled` attribute.

        Args:
            patterns: Mapping an enum representing the states a file could be in
                to a list of regex patterns.
        """
        self._patterns = patterns
        for state, pat_list in patterns.items():
            for pat in pat_list:
                self._compiled.append((re.compile(pat.decode(), re.MULTILINE), state))

    def stream(self, buf: bytes) -> Iterator[tuple[int, int, bytes, S]]:
        text = buf.decode("utf-8", "ignore")
        for rex, state in self._compiled:
            for m in rex.finditer(text):
                yield m.start(), m.end(), m.re.pattern.encode(), state
