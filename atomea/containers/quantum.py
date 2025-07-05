import atomea.typing as adt
from atomea.containers import Container
from atomea.data import Cadence, Data
from atomea.stores import StoreKind


class Quantum(Container):
    """This section encompasses data pertaining to quantum mechanical descriptions."""

    electron_frozen_num = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="5b44b60c-8435-41c4-88d5-cb4a1883b75b",
        description="Total number of frozen electrons",
    )
    """Specifies the total number of electrons considered as frozen in quantum chemical
    calculations.

    Frozen electrons are those that are not included in the active
    space for electronic structure calculations. This approach is used to simplify the
    computational process by reducing the number of electrons that need to be actively
    considered, thereby focusing on those more likely to be involved in chemical
    reactions or significant bonding interactions.
    """

    charge_net = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="6ff82a49-4666-4cbb-978a-409bfa6a511",
        description="Net charge",
    )
    """The net charge of an atomic system is the overall charge determined by the
    balance between positively charged protons and negatively charged electrons.
    """

    multiplicity = Data[adt.DataFrame](
        store_kind=StoreKind.TABLE,
        uuid="8e3eb55a-ed81-46d3-9f34-0ea00fa8c8e4",
        description="Spin multiplicity",
    )
    """The degeneracy or the number of possible spin states associated with a
    particular electronic state of a molecule. The multiplicity is denoted by the
    symbol $2S+1$, where $S$ is the total electron spin angular momentum.
    Here, $S$ can take non-negative half-integer values, such as 0, 1/2, 1, 3/2, and
    so on. The multiplicity is always an integer.
    """

    def __init__(self, parent: object) -> None:
        self.label = "quantum"
        self.cadence = Cadence.ENSEMBLE
        self._parent = parent
        self.electron_frozen_num.bind_to_container(self)
        self.charge_net.bind_to_container(self)
        self.multiplicity.bind_to_container(self)
