from typing import Any

from abc import ABC

import polars as pl


class AtomeaContainer(ABC):
    id: str
    """Label for the container."""

    @classmethod
    def prepare_table(cls, **kwargs: Any) -> pl.DataFrame:
        """Prepare a DataFrame for this container."""
        raise NotADirectoryError
