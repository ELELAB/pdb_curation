# pdb_curation

Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

## Introduction

pdb_curation is a Snakemake pipeline designed to collect protein structures
and ensembles from PDB and process and annotate them in different
ways. 

## Requirements

### Required software

The pipeline uses the following Python packages:

- pandas
- pypdb

### Input structure

The sole input for this pipeline is a `pdbs.csv` file that should be located
in the same directory as the Snakefile (see example included in the repository). 
This file should contain the following comma-separated columns:

|Column name|Expected content|Example|
|------------|----------------|-------|
|pdbs|four-letter PDB ID code. Can be upper or lowercase|1AQH|
|type|denotes if the protein is part of a protein complex (complex) or not (free)|free|

## Output structure

The pipeline considers all the PDB ids in the input csv file and processes them
automatically. One folder per PDB ID is created in the current directory and itse
content varies depending on whether the PDB contains a single protein or a 
protein complex.

Here is an example of the output of the pipeline for the 1AQH PDB, a monomeric
protein:

```
1AQH/
├── 1AQH_original.pdb
└── 1AQH_pdbredo.pdb
pdb_methods.csv
```

- `1AQH/1AQH_original.pdb` is the unprocessed pdb file downloaded from the PDB
database (https://www.rcsb.org/structure/1AQH). 

- `1AQH/1AQH_pdbredo.pdb` is the refined and rebuilt PDB structure from 
the PDB-REDO webserver (e.g. https://pdb-redo.eu/db/1aqh/1aqh_final.pdb)

- `pdb_methods.csv` is a csv table summary file containing the same information 
as the input file and more, namely:

|Column name|Expected content|Example|
|------------|----------------|-------|
|pdb|four-letter PDB ID code. Can be upper or lowercase|1AQH|
|type|denotes if the protein is part of a protein complex (complex) or not (free)|free|
|method|denotes the experimental method used to obtain the structure|X-ray diffraction|

