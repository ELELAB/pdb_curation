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

- pandae
- pypdb

### Input structure

The sole input for this pipeline is a `pdbs.csv` file that should be located
in the same directory as the Snakefile (see example). It includes the following
comma-separated columns:

|Column name|Expected content|Example|
-----------------------------------
pdbs|four-letter PDB ID code. Can be upper or lowercase|1AQH|
...

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

- `1AQH/1AQH_original.pdb` is ..

- `1AQH/1AQH_pdbredo.pdb` is ..

