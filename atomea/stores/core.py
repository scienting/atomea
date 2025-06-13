from abc import ABC

from atomea.stores import DiskFormat


class Store(ABC):
    def __init__(self, disk_format: DiskFormat = DiskFormat.NONE) -> None:
        self.disk_format = disk_format
