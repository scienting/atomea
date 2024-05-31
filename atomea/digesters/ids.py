from typing import Any, Callable

from functools import wraps


def SchemaUUID(uuid_str: str) -> Callable[..., Any]:
    """
    Decorator to assign a specific UUID to a function.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        # Assign the provided UUID to the function
        if not hasattr(func, "__uuid__"):
            func.__uuid__ = uuid_str  # type: ignore
        return func

    return decorator
