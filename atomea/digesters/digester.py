from typing import Any, Generator

from abc import ABC, abstractmethod
from collections.abc import Collection, MutableSequence

from ..schema import Atomea


class Digester(ABC):
    """Digest results into a desired atomea.

    Child classes must implement the digest method to processes all possible
    information.
    """

    @classmethod
    def array_size(
        cls, schema_info: dict[str, Any], *args: Any, **kwargs: Collection[Any]
    ) -> MutableSequence[int]:
        """Return the size of the array to be stored."""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def digest(
        self, atomea: Atomea, *args: Any, **kwargs: Collection[Any]
    ) -> dict[str, Any]:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def digest_step(
        self, atomea: Atomea, *args: Any, **kwargs: Collection[Any]
    ) -> Generator[dict[str, Any], None, None]:
        raise NotImplementedError
