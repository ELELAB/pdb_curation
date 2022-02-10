import pandas as pd
import requests
from urllib import request as ur
from urllib.error import HTTPError
from pypdb import get_all_info
import pypdb

pdb_csv=pd.read_csv("pdbs.csv")

pdb_list_old=pdb_csv['pdb'].to_list()
pdb_list = []
for i in pdb_list_old:
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i+"_pdbredo.pdb")

rule all:
    input:
        pdb_list,
        expand("{pdb}/{pdb}_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
        "pdbs_method.csv"

rule download_pdb:
    output:
        "{pdb}/{pdb}_original.pdb"
    run:
        data=pypdb.get_pdb_file(f'{wildcards.pdb}')
        with open(f"{output}",'w') as f:
            for i in data:
                f.write(i)

rule download_redo:
    output:
        "{pdb}/{pdb_r}_pdbredo.pdb"
    run:
        try:
            ur.urlopen(f"https://pdb-redo.eu/db/{str.lower(wildcards.pdb)}/{str.lower(wildcards.pdb)}_final.pdb")
        except HTTPError:
            pass
        ur.urlretrieve(f"https://pdb-redo.eu/db/{str.lower(wildcards.pdb)}/{str.lower(wildcards.pdb)}_final.pdb", output[0])

rule method_table:
    input:
        "pdbs.csv"
    output:
        "pdbs_method.csv"
    run:
        pdb_csv=pd.read_csv(f"{input}")
        pdb_csv['method'] = pdb_csv.apply(lambda x: get_all_info(x['pdb'])['exptl'][0]['method'], axis=1)
        print(pdb_csv)        
        pdb_csv.to_csv(f"{output}")
             
