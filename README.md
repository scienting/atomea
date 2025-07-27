<h1 align="center">atomea</h1>
<h4 align="center">A hierarchical data model for atomistic systems</h4>
<p align="center">
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

`atomea` requires Python 3.12 or later.
You can currently only install it from the GitHub `main` branch.

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

This software is licensed under the **Prosperity Public License 3.0.0**.
See [LICENSE.md](https://github.com/scienting/atomea/blob/main/LICENSE.md) for full terms, including noncommercial use and 30-day commercial trial conditions.
