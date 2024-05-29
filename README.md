<h1 align="center">atomea</h1>

<h4 align="center">Extensible schema for atomistic simulations and calculations.</h4>

<h4 align="center" style="padding-bottom: 0.5em;"><a href="https://atomea.oasci.org">Documentation</a></h4>

<p align="center">
    <a href="https://gitlab.com/oasci/software/atomea/-/pipelines">
        <img src="https://gitlab.com/oasci/software/atomea/badges/main/pipeline.svg" alt="Build Status ">
    </a>
    <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/atomea">
    <a href="https://codecov.io/gh/oasci/atomea">
        <img src="https://codecov.io/gh/oasci/atomea/branch/main/graph/badge.svg?token=74wLrsOMTD" alt="codecov">
    </a>
    <a href="https://github.com/oasci/atomea/releases">
        <img src="https://img.shields.io/github/v/release/oasci/atomea" alt="GitHub release (latest by date)">
    </a>
    <a href="https://github.com/oasci/atomea/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/oasci/atomea" alt="License">
    </a>
    <a href="https://github.com/oasci/atomea/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/oasci/atomea" alt="GitHub repo size">
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

## Overview

**atomea** is a Python package designed to simplify and standardize the setup and data management for atomistic simulations and calculations.
It leverages the power of Pydantic for schema definition and Jinja2 templates for input file generation, making it easy to automate, document, and prepare input files for various computational chemistry and biology tools.

## Key Features

-   **Extensible Schema Definition:** Easily define schemas for various computational packages using Pydantic.
-   **Automated Input File Preparation:** Generate input files with Jinja2 templates to ensure consistency and reproducibility.
-   **Data Digestion:** Convert raw output files into optimized storage formats (HDF5, Parquet, Zarr) with a consistent interface.
-   **YAML Integration:** Save and load configurations and data in YAML format for easy sharing and reproducibility.

## Installation

You can install **atomea** via pip:

```bash
pip install atomea
```

## Deploying

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings and update [`CHANGELOG.md`](https://gitlab.com/oasci/software/atomea/-/blob/main/CHANGELOG.md).

For example, to bump the `minor` version you would run the following command.

```bash
poetry run bump-my-version bump minor
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## License

This project is released under the Apache-2.0 License as specified in [`LICENSE.md`](https://gitlab.com/oasci/software/atomea/-/blob/main/LICENSE.md).
