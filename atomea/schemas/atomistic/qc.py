from typing import Annotated

from pydantic import BaseModel, Field

from ..io import YamlIO


class QCSchema(BaseModel, YamlIO):
    """This section encompasses data pertaining to quantum mechanical descriptions or
    calculations that are not covered in the [system documentation](./system.md).
    It includes specialized parameters and results specific to quantum mechanics that
    are essential for advanced computational chemistry and physics analysis."""

    electron_frozen_num: Annotated[
        int | None,
        {"cadence": "molecule", "uuid": "5b44b60c-8435-41c4-88d5-cb4a1883b75b"},
    ] = Field(default=None)
    """Specifies the total number of electrons considered as frozen in quantum chemical
    calculations.

    Frozen electrons are those that are not included in the active
    space for electronic structure calculations. This approach is used to simplify the
    computational process by reducing the number of electrons that need to be actively
    considered, thereby focusing on those more likely to be involved in chemical
    reactions or significant bonding interactions.

    **Cadence:** `molecule`

    **UUID:** `5b44b60c-8435-41c4-88d5-cb4a1883b75b`
    """
