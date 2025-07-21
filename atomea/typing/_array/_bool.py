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

Bool: TypeAlias = npt.NDArray[np.bool_]

OptionalBool: TypeAlias = Bool | None
