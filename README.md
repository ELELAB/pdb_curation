Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark

Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

pdb_curation

pdb_curation is a pipeline Snakemake based designed to collect protein structures and ensembles from pdb or pdb_redo databases and parse them for further analysis. The structures are downloaded and sorted in different folders, accordingly to the database, along with a file that collects the methods applied to obtain the structures. The structure sequences are then aligned with ClustalOmega against the corrispettive UNIPROT sequence in order to retrive all the missing atoms and residues. A following Modeller step will be applied to add all the missing part in the pdb files. It includes Python libraries and ClustalOmega and Modeller softwares.

