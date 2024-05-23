from pydantic import BaseModel, Field


class IdentificationSchema(BaseModel):
    """Information that assists in uniquely identifying this data."""

    sha256_coordinates: str | None = Field(default=None)
    """A sha256 hash based only on coordinates."""
