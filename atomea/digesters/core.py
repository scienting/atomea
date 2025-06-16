from typing import Any

from abc import ABC, abstractmethod

from loguru import logger

from atomea.containers import Project


class Digester(ABC):
    """
    Base class for all digesters.
    """

    @classmethod
    @abstractmethod
    def prepare(cls, *args: Any, **kwargs: Any) -> dict[str, Any]:
        """
        Load files, open handles, build FSM state, etc.
        Return a dict representing the full parsing context.
        """
        ...

    @classmethod
    def checks(cls) -> None:
        """
        Optional pre‐flight checks (e.g. dependencies).
        Override if needed.
        """
        pass

    @classmethod
    @abstractmethod
    def extract(
        cls, prj: Project, id_ens: str, id_run: str, ctx: dict[str, Any]
    ) -> Project:
        """Extract and parse all possible information given the context.
        This method is responsible for calling all other methods.

        Args:
            prj: Destination project to put extracted information.
            id_ens: ID to store any ensemble data under.
            context: Information needed for the digestion process.

        Returns:
            Project after digestion.
        """
        ...

    @classmethod
    def run(
        cls,
        prj: Project,
        id_ens: str,
        id_run: str,
        digest_args: tuple[Any] | None = None,
        digest_kwargs: dict[str, Any] | None = None,
    ) -> Project:
        """
        Run the digestion process from start to finish for a single ensemble.

        Args:
            prj: Project to store all digested data to.
            id_ens: ID of this ensemble. This function will create the ensemble
                if it does not exist in `prj`.
            digest_args: Arguments for preparing the context needed for digestion.
            digest_kwargs: Keyword arguments for preparing the context needed for
                the digestion.

        Returns:
            Project after digestion.
        """
        logger.info("Running digestion")
        if digest_args is None:
            digest_args = tuple()
        if digest_kwargs is None:
            digest_kwargs = dict()
        cls.checks()
        ctx = cls.prepare(*digest_args, **digest_kwargs)

        _ = prj.get_ensemble(id_ens) or prj.add_ensemble(id_ens)

        project = cls.extract(prj, id_ens, id_run, ctx)
        for store_kind, store in project._stores.items():
            store.dump()

        return project
