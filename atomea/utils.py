import importlib


def get_obj_from_string(import_string):
    """Retrieves an object based on an import string and object name.

    Args:
        import_string: The import string, starting from the root module, containing
            the desired object. This function would be
            `"atomea.utils.get_obj_from_string"`.

    Returns:
        The object identified by the import string.
    """
    # Try to import as object inside a module.
    try:
        module_name, obj_name = import_string.rsplit(".", 1)
        module = importlib.import_module(module_name)
        obj = getattr(module, obj_name)
    # Fallback to trying the import as a module.
    except AttributeError:
        obj = importlib.import_module(import_string)
    return obj
