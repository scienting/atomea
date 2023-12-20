#!/usr/bin/env python3

import os

import yaml

SRC_DIR = "atomea/schemas"
DEST_DIR = "docs/schemas"


# Collect paths of all YAML files
def get_files(path, expression, recursive=True):
    r"""Returns paths to all files in a given directory that matches a provided
    expression in the file name.

    Commonly used to find all files of a certain type, e.g., output or xyz
    files.

    Parameters
    ----------
    path : :obj:`str`
        Specifies the directory to search.
    expression : :obj:`str`
        Expression to be tested against all file names in ``path``.
    recursive : :obj:`bool`, default: ``True``
        Recursively find all files in all subdirectories.

    Returns
    -------
    :obj:`list` [:obj:`str`]
        All absolute paths to files matching the provided expression.
    """
    # pylint: disable=redefined-outer-name
    if path[-1] != "/":
        path += "/"
    if recursive:
        all_files = []
        for dirpath, _, filenames in os.walk(path):
            index = 0
            while index < len(filenames):
                if dirpath[-1] != "/":
                    dirpath += "/"
                filenames[index] = dirpath + filenames[index]
                index += 1
            all_files.extend(filenames)
        files = []
        for f in all_files:
            if expression in f:
                files.append(f)
    else:
        files = []
        for f in os.listdir(path):
            filename = os.path.basename(f)
            if expression in filename:
                files.append(path + f)
    return files


yaml_file_paths = get_files(SRC_DIR, ".yml")
yaml_file_paths = [os.path.abspath(path) for path in yaml_file_paths]


# Load all YAML files.
yaml_files = {}
for file_path in yaml_file_paths:
    file_name = os.path.splitext(os.path.basename(file_path))[0]

    with open(file_path, encoding="utf-8") as f:
        yaml_files[file_name] = yaml.safe_load(f)

for file_name, defs in yaml_files.items():
    file_path = os.path.join(DEST_DIR, file_name + ".md")

    with open(file_path, mode="w", encoding="utf-8") as f:
        f.write("# " + file_name.capitalize() + "\n")

        for defs_cat_name, defs_dicts in defs.items():
            f.write("\n## `" + defs_cat_name + "`\n")

            f.write("\n" + defs_dicts["description"])

            f.write("\n**Dimensions:** " + str(defs_dicts["ndim"]) + "\n")
            f.write("\n**Data type:** `" + defs_dicts["dtype"] + "`\n")
            if defs_dicts["units"] is not None:
                f.write("\n**Units:** " + defs_dicts["units"] + "\n")
