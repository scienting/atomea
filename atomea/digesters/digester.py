from typing import Any, Generator

from abc import ABC, abstractmethod
from collections.abc import Collection

from ..schema import Atomea


class Digester(ABC):
    """Digest results into a desired atomea.

    Child classes must implement the digest method to processes all possible
    information.
    """

    @classmethod
    @abstractmethod
    def digest(
        self, atomea: Atomea, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def digestStep(
        self, atomea: Atomea, *args: Any, **kwargs: Collection[Any]
    ) -> Generator[dict[str, Any], None, None]:
        raise NotImplementedError
