from typing import Iterator

from dataclasses import dataclass


@dataclass
class FileChunk:
    """A chunk of file data with position info"""

    data: bytes
    start_pos: int
    end_pos: int
    is_last: bool = False


class ChunkReader:
    """Manages reading file in chunks"""

    def __init__(self, file_path: str, chunk_size: int = 65536, overlap: int = 1024):
        self.file_path = file_path
        self.chunk_size = chunk_size
        self.overlap = overlap

    def iter_chunks(self) -> Iterator[FileChunk]:
        """Iterate over file chunks with overlap"""
        ...
