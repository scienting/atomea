import numpy as np
import numpy.typing as npt

from atomea.schemas.field import Cadence, FieldMeta, SchemaField, StoreKind


class Quantum:
    """This section encompasses data pertaining to quantum mechanical descriptions or
    calculations that are not covered in the [system documentation](./microstates.md).
    It includes specialized parameters and results specific to quantum mechanics that
    are essential for advanced computational chemistry and physics analysis."""

    electron_frozen_num: npt.NDArray[np.uint8] | None = SchemaField[
        npt.NDArray[np.uint8]
    ](
        dtype=np.uint8,
        meta=FieldMeta(
            uuid="5b44b60c-8435-41c4-88d5-cb4a1883b75b",
            cadence=Cadence.ENSEMBLE,
            store=StoreKind.TABLE,
            description="Total number of frozen electrons",
        ),
        default=None,
    )
    """Specifies the total number of electrons considered as frozen in quantum chemical
    calculations.

    Frozen electrons are those that are not included in the active
    space for electronic structure calculations. This approach is used to simplify the
    computational process by reducing the number of electrons that need to be actively
    considered, thereby focusing on those more likely to be involved in chemical
    reactions or significant bonding interactions.

    **Cadence:** `ENSEMBLE`

    **UUID:** `5b44b60c-8435-41c4-88d5-cb4a1883b75b`
    """
