# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import TypeAlias

import numpy as np
import numpy.typing as npt
import polars as pl

PassableData: TypeAlias = npt.NDArray[np.generic] | pl.DataFrame
"""Core data types that we stick to passing around functions. All functions working
with data should accept and/or return data of this type.

All stores should return one of these types.
"""
OptionalPassableData: TypeAlias = npt.NDArray[np.generic] | pl.DataFrame | None
