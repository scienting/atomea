from abc import ABC, abstractmethod
from typing import Any
import polars as pl


class AtomeaContainer(ABC):
    @abstractmethod
    @classmethod
    def prepare_table(cls, **kwargs: Any) -> pl.DataFrame:
        """Prepare a DataFrame for this container."""
        raise NotADirectoryError
