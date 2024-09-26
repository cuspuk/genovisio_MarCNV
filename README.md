# Genovisio MarCNV

[![Python version](https://img.shields.io/badge/python-3.12+-green.svg)](https://www.python.org/downloads/)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

Automatic evaluation of ACMG criteria for CNV based on .json file with data from MongoDB.

## Installation

In python3.12 you can simply use pip to install marCNV:

```bash
pip install git+https://github.com/cuspuk/genovisio_MarCNV.git
```

Without python3.12, you can install marcnv using mamba:

```bash
mamba env create -f conda_marcnv.yaml
```

This gives you 2 entrypoints:

- `marcnv-batch` - running marCNV in batch-processing mode
- `marcnv-single` - running marCNV in single CNV processing mode

## Running

To run marCNV, running instance of mongo database is required. Mongo URI and database name can be supplied to the entrypoint commands, see `--help`. Default MongoDB URI is `mongodb://localhost:27017/` and the database name 'genovisio'.

To run marCNV, call one of entrypoint commands (if installed using conda, activate it first), like:

```sh
marcnv-single chr16:34289161-34490212/loss
```

### Options

Adding flag `-j` can output json file with all needed data for classification.
Flag `-s` skips the automatic classification if you need only the json output.

```bash
marcnv-single chr16:34289161-34490212/loss -j test.json -s
```

### Batch annotation

Alternatively use the batch script if you have more CNVs needed to classify (accepts CNVs in .bed or .tsv file, outputs .tsv file with all annotated
information).

```bash
marcnv-batch cnvs.bed cnvs_annotated.tsv
```

## Development

There are some changes whether you develop on machines with Python 3.12 and on machines without Python 3.12.
If you cannot install python3.12 and are limited to conda environments, skip to the next guide.

### Development with Python 3.12

Install poetry using pipx:

```sh
pipx install poetry
```

Now in the cloned repository, install the package:

```sh
poetry install
```

Activate the virtual environment where dependencies are installed:

```sh
poetry shell
```

All dependencies are now installed and you can run any marcnv code using:

- entrypoint commands `marcnv-single {ARGS}` and `marcnv-batch {ARGS}`
- `python {python_script}` like `python marcnv/classify_cnv.py {ARGS}`

### Development without python 3.12

Install custom conda environment with python3.12.

Then, install poetry (python packaging management library):

```sh
curl -sSL https://install.python-poetry.org | python3 -
```

Now in the cloned repository, install the package:

```sh
poetry install
```

All dependencies are now installed and you can run any marcnv code using:

- entrypoint commands `marcnv-single {ARGS}` and `marcnv-batch {ARGS}`
- `python {python_script}` like `python marcnv/classify_cnv.py {ARGS}`

### Adding/removing dependencies

When adding or removing dependencies, you need to define them in `pyproject.toml`.

Then, locked versions need to be redefined:

```sh
poetry lock
```

Install again:

```sh
poetry install
```

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.
This means that there the PR title must be like `feat: XY` or `fix: XY` and so on, and the PR should contain at least one commit named like this.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history
