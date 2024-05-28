from pydantic import BaseModel, Field

from ..io import IOBase
from .cfg.main import WestpaConfigConfig
from .environment import WestpaEnv


class WestpaConfig(BaseModel, IOBase):
    """The root configuration for `west.cfg` files."""

    west: WestpaConfigConfig = Field(default_factory=WestpaConfigConfig)
    env: WestpaEnv = Field(default_factory=WestpaEnv)
