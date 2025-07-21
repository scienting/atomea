# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

from abc import ABC

import polars as pl


class Container(ABC):
    label: str
    """Label for the container."""

    @classmethod
    def prepare_table(cls, **kwargs: Any) -> pl.DataFrame:
        """Prepare a DataFrame for this container."""
        raise NotImplementedError
