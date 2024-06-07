from pydantic import BaseModel, Field

from .....render import Render
from .data import DataConfig
from .executable import ExecutableConfig
from .propagation import PropagationConfig
from .system import SystemConfig


class WestpaConfigConfig(BaseModel, Render):
    """
    Configuration class for wrapping WESTPA configuration options. These options
    are scattered across a variety of sources shown below. We do our best to aggregate
    and explain each parameter.

    -   [GitHub wiki](https://github.com/westpa/westpa/wiki/Configuration-File),
    -   [readthedocs](https://westpa.readthedocs.io/en/latest/users_guide/west/setup.html#configuration-file),

    Attributes:
        system (SystemConfig): Configuration options for the system.
        propagation (PropagationConfig): Configuration options for the propagation.
        data (DataConfig): Configuration options for the data.
        executable (ExecutableConfig): Configuration options for the executable.
    """

    system: SystemConfig = Field(default_factory=SystemConfig)
    propagation: PropagationConfig = Field(default_factory=PropagationConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    executable: ExecutableConfig = Field(default_factory=ExecutableConfig)
