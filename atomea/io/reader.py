# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

from abc import ABC, abstractmethod

from loguru import logger

from atomea.containers import Project


class Reader(ABC):
    """
    Base class for all readers.
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
        Optional pre‐flight checks (e.g. dependencies).
        Override if needed.
        """

    @classmethod
    @abstractmethod
    def extract(
        cls, prj: Project, ens_id: str, run_id: str, ctx: dict[str, Any]
    ) -> Project:
        """Extract and parse all possible information given the context.
        This method is responsible for calling all other methods.

        Args:
            prj: Destination project to put extracted information.
            ens_id: ID to store any ensemble data under.
            context: Information needed for the reader.

        Returns:
            Project after extracting all information.
        """

    @classmethod
    def run(
        cls,
        prj: Project,
        ens_id: str,
        run_id: str,
        reader_args: tuple[Any] | None = None,
        reader_kwargs: dict[str, Any] | None = None,
    ) -> Project:
        """
        Read and extract data from start to finish for a single ensemble.

        Args:
            prj: Project to store all extracted data to.
            ens_id: ID of this ensemble. This function will create the ensemble
                if it does not exist in `prj`.
            reader_args: Arguments for preparing the context needed for the reader.
            reader_kwargs: Keyword arguments for preparing the context needed for
                the reader.

        Returns:
            Project after extracting data.
        """
        logger.info("Running the reader")
        if reader_args is None:
            reader_args = tuple()
        if reader_kwargs is None:
            reader_kwargs = dict()
        cls.checks()
        ctx = cls.prepare(*reader_args, **reader_kwargs)

        _ = prj.get_ensemble(ens_id) or prj.add_ensemble(ens_id)

        project = cls.extract(prj, ens_id, run_id, ctx)
        for store_kind, store in project._stores.items():
            store.flush()

        return project
