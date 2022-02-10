import pandas as pd
import requests
from urllib import request as ur
from urllib.error import HTTPError
from pypdb import get_all_info

pdb_csv=pd.read_csv("pdbs.csv")

pdb_list_old=pdb_csv['pdb'].to_list()
pdb_list = []
for i in pdb_list_old:
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i+"_final.pdb")

rule all:
    input:
        pdb_list,
        expand("{pdb}/{pdb}.pdb", pdb=pdb_csv['pdb'].str.upper())
        
rule download_pdb:
    output:
        "{pdb}/{pdb}.pdb"
    shell:
        """
        wget https://files.rcsb.org/view/{wildcards.pdb}.pdb -P {wildcards.pdb}
        """

rule download_redo:
    output:
        "{pdb}/{pdb_r}_final.pdb"
    run:
        for i in pdb_list:
            i = i.lower()
            shell("wget https://pdb-redo.eu/db/"+i+" -P {wildcards.pdb}") 
