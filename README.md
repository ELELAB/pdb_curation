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
- urllib

### Input structure

The sole input for this pipeline is a `pdbs.csv` file that should be located
in the same directory as the Snakefile (see example included in the repository). 
This file should contain the following comma-separated columns:

|Column name|Expected content|Example|
|------------|----------------|-------|
|pdbs|four-letter PDB ID code. Can be upper or lowercase|1AQH|
|type|denote if the protein is in complex with a ligand (complex) or not (free)|free|


## Output structure

The pipeline considers all the PDB id in the input csv file and processes them
automatically. One folder per PDB ID is created in the current directory whose
content varies depending on whether the PDB contains a single protein or a 
protein complex.

Here is an example of the output of the pipeline for the 1AQH PDB, a monomeric
protein:

```
1AQH/
├── 1AQH_original.pdb
└── 1AQH_pdbredo.pdb
```

- `1AQH/1AQH_original.pdb` is the path containing the unprocessed pdb file
   downloaded from PDB database (https://www.rcsb.org/structure/1AQH). 

- `1AQH/1AQH_pdbredo.pdb` is the path containing the pdb file optimized through
   the PDB-REDO procedure (https://pdb-redo.eu/db/1aqh/1aqh_final.pdb).

For each pdb entry contained in the 'pdb.csv' file, the pipeline extracts the
information about 'method' (how the structure has been solved)  and 'type' (if the
structure is in complex with a ligand or not) and sorts them in a csv file (pdb_methods.csv) located in the same folder of the Snakefile and input files.

The output file will contain the following comma-spearated columns:

|Column name|Expected content|Example|
|------------|----------------|-------|
|pdb|four-letter PDB ID code. Can be upper or lowercase|1AQH|
|type|denote if the protein is in complex with a ligand (complex) or not (free)|free|
|method|denote the type of experiment performed to obtain the structure|X-ray diffraction|



