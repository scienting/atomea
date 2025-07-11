from typing import TypeAlias

import numpy as np
import numpy.typing as npt

Bool: TypeAlias = npt.NDArray[np.bool_]

OptionalBool: TypeAlias = Bool | None
