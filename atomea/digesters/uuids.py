from typing import Any, Callable

from functools import wraps


# pylint: disable-next=invalid-name
def SchemaUUID(uuid_str: str) -> Callable[..., Any]:
    """
    Decorator to assign a specific Universally Unique Identifier (UUID) to a method,
    particularly a digester method for extracting data from frame inputs.

    The UUID is stored in the private attribute `__uuid__` of the method, which can
    then be used to map the method to specific fields in a schema
    (e.g., [`MoleculeSchema`][schemas.atomistic.MoleculeSchema]).

    This is particularly useful when dealing with data extraction and transformation
    processes in computational chemistry and biology. By using
    [`@SchemaUUID`][digesters.ids.SchemaUUID], each
    method that processes a part of the input data can be easily identified and called
    based on its UUID. This allows for a flexible and dynamic way to handle various
    data processing tasks, ensuring that each piece of data is processed by the
    appropriate method.

    Args:
        uuid_str: The UUID to be assigned to the method.

    Returns:
        The decorated function with the assigned UUID.

    Example:
        ```python
        @SchemaUUID("81c7cec9-beec-4126-b6d8-91bee28951d6")
        def coordinates(atoms: mda.AtomGroup) -> npt.NDArray[np.float64]:
            # Method implementation
        ```

    Notes:
        - The decorator ensures that the UUID is only assigned once to avoid
          accidental overwriting.
        - The UUID can be retrieved using the `__uuid__` attribute of the decorated method.
        - This decorator is particularly useful for organizing and identifying methods
          that process specific parts of input data. When combined with a mapping
          function like `get_uuid_map`, it allows for dynamic method calling based on
          the data being processed.
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
