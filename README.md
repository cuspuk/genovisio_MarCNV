# Genovisio MarCNV

Automatic evaluation of ACMG criteria for CNV based on .json file with data from MongoDB.

## Installation

Create conda repository with all necessary packages:

```bash
mamba env create -n genovisio_marcnv -f conda_env.yaml
```

### Run mongo

MongoDB should be running before you run the MarCNV:

```bash
conda activate genovisio_marcnv
mongod --dbpath <dbpath_dirname>
```

Where `<dbpath_dirname>` should be an existing directory filled with databases using the
repository [genovisio_sources](ttps://github.com/marcelTBI/genovisio_sources).

## Running the annotation

Activate the conda environment:

```bash
conda activate genovisio_marcnv
```

### Annotation of a single CNV

Run classification script to classify your target CNV (in this case `chr16:34289161-34490212/loss`):

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
