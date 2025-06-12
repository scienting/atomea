from typing import Any, Collection

import inspect
from abc import ABC, abstractmethod

from loguru import logger

from atomea.containers import Project
from atomea.stores.arrays import ArrayStore
from atomea.stores.tables import TableStore


class Digester(ABC):
    """
    Base class for all digesters.  Subclasses must define:
      - prepare_inputs(...) → context
      - any number of `parse_<thing>(context, project, ensemble_id)` methods
    """

    @classmethod
    @abstractmethod
    def prepare_inputs(cls, *args: Any, **kwargs: Collection[Any]) -> dict[str, Any]:
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
    def digest(
        cls,
        arrays: ArrayStore,
        tables: TableStore,
        *prepare_args: Any,
        project: Project | None = None,
        ensemble_id: str = "default",
        **prepare_kwargs: Any,
    ) -> Project:
        """
        Top‐level entrypoint: builds or reuses a Project;
        creates or gets an Ensemble by `ensemble_id`;
        calls prepare_inputs, then runs every parse_* method.
        """
        cls.checks()
        ctx = cls.prepare_inputs(*prepare_args, **prepare_kwargs)

        if project is None:
            project = Project(arrays, tables)

        # ensure the requested ensemble exists
        ens = project.get_ensemble(ensemble_id) or project.add_ensemble(ensemble_id)

        # discover all parse_ methods
        for name, method in inspect.getmembers(cls, predicate=inspect.isfunction):
            if not name.startswith("parse_"):
                continue
            logger.debug(f"[{ensemble_id}] running parser `{name}`")
            try:
                # each parse_ method takes (context, project, ensemble_id)
                method(ctx, project, ensemble_id)
            except Exception:
                logger.exception(f"[{ensemble_id}] error in parser `{name}`")

        return project
