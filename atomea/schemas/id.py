from typing import Annotated

from pydantic import BaseModel, Field


class IdentificationSchema(BaseModel):
    """Information that assists in uniquely identifying this data."""

    sha256_coordinates: Annotated[
        str | None,
        {"cadence": "ensemble", "uuid": "a6a84df4-f182-4e4f-bd52-168004442ec2"},
    ] = Field(default=None)
    """A sha256 hash based only on coordinates.

    **Cadence:** `ensemble`

    **UUID:** `a6a84df4-f182-4e4f-bd52-168004442ec2`
    """
