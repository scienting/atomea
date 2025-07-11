from typing import Any

from abc import ABC, abstractmethod

from loguru import logger

from atomea.containers import Project


class Writer(ABC):
    """
    Base class for all writers.
    """

    @classmethod
    @abstractmethod
    def prepare(cls, *args: Any, **kwargs: Any) -> dict[str, Any]:
        """
        Load files, open handles, build FSM state, etc.
        Return a dict representing the full parsing context.
        """

    @classmethod
    def checks(cls) -> None:
        """
        Optional pre-flight checks (e.g. dependencies).
        Override if needed.
        """

    @classmethod
    def run(
        cls,
        prj: Project,
        ens_id: str,
        run_id: str,
        writer_args: tuple[Any] | None = None,
        writer_kwargs: dict[str, Any] | None = None,
    ) -> Project:
        """
        Write data for a single ensemble run.

        Args:
            prj: Project to store all extracted data to.
            ens_id: ID of this ensemble. This function will create the ensemble
                if it does not exist in `prj`.
            writer_args: Arguments for preparing the context needed for the writer.
            writer_kwargs: Keyword arguments for preparing the context needed for
                the writer.
        """
        logger.info("Running the writer")
        if writer_args is None:
            writer_args = tuple()
        if writer_kwargs is None:
            writer_kwargs = dict()
        cls.checks()
        ctx = cls.prepare(*writer_args, **writer_kwargs)

        _ = prj.get_ensemble(ens_id) or prj.add_ensemble(ens_id)

        # TODO: implement
