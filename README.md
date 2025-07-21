<h1 align="center">atomea</h1>
<h4 align="center">A hierarchical data model for atomistic calculations</h4>
<p align="center">
    <a href="https://github.com/scienting/atomea/actions/workflows/tests.yml">
        <img src="https://github.com/scienting/atomea/actions/workflows/tests.yml/badge.svg" alt="Build Status ">
    </a>
    <a href="https://github.com/scienting/atomea/actions/workflows/docs.yml">
        <img src="https://github.com/scienting/atomea/actions/workflows/docs.yml/badge.svg" alt="Build Status ">
    </a>
    <a href="https://badge.fury.io/py/atomea">
        <img src="https://badge.fury.io/py/atomea.svg" alt="PyPI version" height="18">
    </a>
    <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/atomea">
    <a href="https://codecov.io/gh/scienting/atomea">
        <img src="https://codecov.io/gh/scienting/atomea/branch/main/graph/badge.svg" alt="codecov">
    </a>
    <a href="https://github.com/scienting/atomea/releases">
        <img src="https://img.shields.io/github/v/release/scienting/atomea" alt="GitHub release (latest by date)">
    </a>
    <a href="https://github.com/scienting/atomea/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/scienting/atomea" alt="License">
    </a>
    <a href="https://github.com/scienting/atomea/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/scienting/atomea" alt="GitHub repo size">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
    <a href="https://github.com/PyCQA/pylint" target="_blank">
        <img src="https://img.shields.io/badge/linting-pylint-yellowgreen" alt="Black style">
    </a>
    <a href="https://github.com/astral-sh/ruff" target="_blank">
        <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Black style">
    </a>
</p>
<h4 align="center"><a href="https://atomea.scient.ing">Documentation</a></h4>

## Overview

Atomea is a Python package for atomistic data and calculations.

## Installation

`atomea` requires Python 3.10 or later.
You can install it using pip, the Python package installer,

```bash
pip install atomea
```

or directly from the GitHub `main` branch.

```bash
pip install git+https://github.com/scienting/atomea.git
```

## Development

We use [pixi](https://pixi.sh/latest/) to manage Python environments and simplify the developer workflow.
Once you have [pixi](https://pixi.sh/latest/) installed, move into `atomea` directory (e.g., `cd atomea`) and install the  environment using the command

```bash
pixi install
```

Now you can activate the new virtual environment using

```sh
pixi shell
```

## License

This project is released under the Apache-2.0 License as specified in [`LICENSE.md`](https://gitlab.com/scienting/software/atomea/-/blob/main/LICENSE.md).

<!-- REFERENCES -->

[zarr]: https://zarr.dev/
