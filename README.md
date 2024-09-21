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
mamba env create -f marcnv.yaml
```

This gives you 2 entrypoints:

- `marcnv-batch` - running marCNV in batch-processing mode
- `marcnv-single` - running marCNV in single CNV processing mode

## Running

To run marCNV, running instance of mongo database is required. Mongo URI and database name can be supplied to the entrypoint commands, see `--help`. Default MongoDB URI is `mongodb://localhost:27017/` and the database name 'genovisio'.

To run marCNV, call one of entrypoint commands (if installed using conda, activate it first).

### Annotation of a single CNV

To classify your target CNV (in this case `chr16:34289161-34490212/loss`) run:

```bash
python classify_cnv.py chr16:34289161-34490212/loss
```

The same script with flag `-j` can output json file with all needed data for classification.
Flag `-s` skips the automatic classification if you need only the json output

```bash
python classify_cnv.py chr16:34289161-34490212/loss -j test.json -s
```

### Batch annotation

Alternatively use the batch script if you have more CNVs needed to classify (accepts CNVs in .bed or .tsv file, outputs .tsv file with all annotated
information).

```bash
python classify_batch.py cnvs.bed cnvs_annotated.tsv
```
