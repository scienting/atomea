from pydantic import BaseModel, Field

from atomea.schemas.io import YamlIO
from atomea.schemas.workflow.westpa import WestpaEnv
from atomea.schemas.workflow.westpa.cfg.core import WestpaConfigConfig


class WestpaConfig(BaseModel, YamlIO):
    """The root configuration for `west.cfg` files."""

    west: WestpaConfigConfig = Field(default_factory=WestpaConfigConfig)
    env: WestpaEnv = Field(default_factory=WestpaEnv)
