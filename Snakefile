configfile: "Config.yaml"
import pandas as pd
import requests
from urllib import request as ur
from urllib.error import HTTPError
from pypdb import get_all_info
import pypdb
import os
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

pdb_csv=pd.read_csv("pdbs.csv")

### uniprot_chain function to extract from each pdb entry the chain and the corrispective uniprot ID ### 

def uniprot_chain(pdb_entry):
    import requests as rq
    import json
    from urllib import request as ur
    url = 'https://data.rcsb.org/graphql'
    headers = { "content-type": "application/graphql" }

    pdb = pdb_entry

    req = f'''{{ entry(entry_id : "{pdb}") {{
        polymer_entities {{
          rcsb_id
          rcsb_polymer_entity_container_identifiers {{
            entity_id
            reference_sequence_identifiers {{
              database_accession
              database_name
            }}
          }}
        entity_poly {{
          rcsb_entity_polymer_type
         }}

        polymer_entity_instances {{
          rcsb_polymer_entity_instance_container_identifiers {{
            auth_asym_id
            }}
          }}
        }}
      }}
    }}'''

    final_correspondence = {}
    response = rq.post(url, headers=headers, data=req)
    entry = json.loads(response.text)
    for entity in entry['data']['entry']['polymer_entities']:
        if entity['entity_poly']['rcsb_entity_polymer_type'] == 'Protein':
            uniprot=[]
            for i in entity['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers']:
                uniprot.append(i['database_accession'])
                for instance in entity['polymer_entity_instances']:
                    final_correspondence[instance['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']] = uniprot
    return final_correspondence

#### list of pdb for which the pdb_redo is available #### 
pdb_list = []
ID_entry_redo=[]
for i in pdb_csv['pdb'].str.lower().to_list():
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i.upper()+"_pdbredo.pdb")
    ID_entry_redo.append(i.upper())

##### pdb_method file #####
pdb_csv['method'] = pdb_csv.apply(lambda x: get_all_info(x['pdb'])['exptl'][0]['method'], axis=1)
pdb_csv.to_csv("pdb_method")
list_chain_ID=[]
for i in pdb_csv['pdb'].to_list():
    dictionary=uniprot_chain(str(i))
    ID_a=", ".join(f"{k} {v}" for k,v in dictionary.items())
    list_chain_ID.append(ID_a)
for i in list_chain_ID:
    pdb_csv["chain_uniprot_ID"]=list_chain_ID
pdb_csv.to_csv("pdb_method")

###dataframes containing the chain ID for each pdb file along with the corresponding uniprot ID:
df=pd.DataFrame()
pdb_id=[]
for entry in pdb_csv['pdb']:
    uniprot_identifiers=uniprot_chain(entry)
    for uniprot_id_list in uniprot_identifiers.values():
        for uniprot_id in uniprot_id_list:
            pdb_id.append(entry.upper())
df['pdb']=pdb_id 
chain_id_list=[]
uniprot_id_final_list=[]
for i in pdb_csv['pdb']:
    uniprot_identifiers=uniprot_chain(i)
    for chain_id,uniprot_id_list in uniprot_identifiers.items():
       for i,g in zip(chain_id,uniprot_id_list):
           for uniprot_id in uniprot_id_list:
               chain_id_list.append(chain_id)
               uniprot_id_final_list.append(uniprot_id)
identifiers=list(zip(chain_id_list,uniprot_id_final_list))       
df1 =pd.DataFrame(identifiers, columns= ['chain_id',"uniprot_id"])     
df_final=pd.concat([df, df1], axis=1, join="inner")       

### dataframe for uniprot_sequence rule:
df_uniprot=df_final.drop_duplicates(subset=['uniprot_id']+['pdb'], keep="last")


### dataframe for split_chain_original rule:
df_split_chain=df_final.drop_duplicates(subset=['chain_id']+['pdb'], keep="last")

###dataframes containing the chain ID for each pdb file for which the pdb redo is available along with the corresponding uniprot ID:

pdb_list = []
ID_entry_redo=[]
for i in pdb_csv['pdb'].str.lower().to_list():
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i.upper()+"_pdbredo.pdb")
    ID_entry_redo.append(i.upper())

df_final_redo=df_final[df_final.pdb.isin(ID_entry_redo)]

### dataframe for split_chain_redo rule:
df_split_chain_redo=df_final_redo.drop_duplicates(subset=['chain_id']+['pdb'], keep="last")




rule all:
    input:
        expand("{pdb}/{pdb}_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb_redo}/{pdb_redo}_redo.pdb", pdb_redo=ID_entry_redo),
        expand("{pdb_align}/alignment/original/{pdb_align}_{chain_id}_{uniprot_seq}.clu", zip, pdb_align=df_split_chain['pdb'], chain_id=df_split_chain['chain_id'], uniprot_seq=df_split_chain['uniprot_id']),
        expand("{pdb_align}/alignment/redo/{pdb_align}_{chain_id}_{uniprot_seq}.clu", zip, pdb_align=df_split_chain_redo['pdb'], chain_id=df_split_chain_redo['chain_id'], uniprot_seq=df_split_chain['uniprot_id']),
        expand("{pdb}/pdb4amber/original/{pdb}_amber_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb_redo}/pdb4amber/redo/{pdb_redo}_amber_redo.pdb",  pdb_redo=ID_entry_redo),


rule download_pdb:
    output:
        "{pdb}/{pdb}_original.pdb",
    run:
        with open(f"{output}",'w') as f:
            f.write(pypdb.get_pdb_file(f'{wildcards.pdb}'))

rule download_redo:
    output:
        "{pdb_redo}/{pdb_redo}_redo.pdb"
    run:
        ur.urlretrieve(f"https://pdb-redo.eu/db/{str.lower(wildcards.pdb_redo)}/{str.lower(wildcards.pdb_redo)}_final.pdb", output[0])



rule uniprot_sequence:
    output:
        "{pdb_uniprot}/uniprot_sequences/{uniprot_seq}.fasta"
    run:
        ur.urlretrieve(f'https://www.uniprot.org/uniprot/{wildcards.uniprot_seq}.fasta', output[0])

rule split_chain_original:
    input:
       "{pdb_split}/{pdb_split}_original.pdb",
    output:
       "{pdb_split}/split_chain/original/{pdb_spli}_{chain_id}_{uniprot_seq}.pdb"
    run:
       structure=PDBParser().get_structure('protein', input[0])
       io=PDBIO()
       for i in structure.get_chains():
           io.set_structure(i)
           io.save(output[0])

rule split_chain_redo:
    input:
       "{pdb_split}/{pdb_split}_redo.pdb",
    output:
       "{pdb_split}/split_chain/redo/{pdb_split}_{chain_id}_{uniprot_seq}.pdb"
    run:
       structure=PDBParser().get_structure('protein', input[0])
       io=PDBIO()
       for i in structure.get_chains():
           io.set_structure(i)
           io.save(output[0])


rule pdb2fasta_original:
     input:
        "{pdb_split}/split_chain/original/{pdb_spli}_{chain_id}_{uniprot_seq}.pdb"
     output:
        "{pdb_split}/split_chain/original/{pdb_spli}_{chain_id}_{uniprot_seq}.fasta"
     shell:
        """
        {config[pdb2fasta_bin]} {input} > {output}
        """



rule pdb2fasta_redo:
     input:
        "{pdb_split}/split_chain/redo/{pdb_spli}_{chain_id}_{uniprot_seq}.pdb"
     output:
        "{pdb_split}/split_chain/redo/{pdb_spli}_{chain_id}_{uniprot_seq}.fasta"
     shell:
        """
        {config[pdb2fasta_bin]} {input} > {output}
        """

rule alignement:
     input:
        "{pdb_align}/split_chain/original/{pdb_align}_{chain_id}_{uniprot_seq}.fasta",
        "{pdb_align}/uniprot_sequences/{uniprot_seq}.fasta"
     output:
        "{pdb_align}/alignment/original/{pdb_align}_{chain_id}_{uniprot_seq}.clu"
     run:
        if str(wildcards.uniprot_seq) in str(input[1]):
            shell("cat {wildcards.pdb_align}/split_chain/original/{wildcards.pdb_align}_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta {wildcards.pdb_align}/uniprot_sequences/{wildcards.uniprot_seq}.fasta > {wildcards.pdb_align}/alignment/original/input_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta")
            shell("{config[clustlo_bin]} -i {wildcards.pdb_align}/alignment/original/input_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta -o {output} --outfmt=clustal --resno")

rule alignement_redo:
     input:
        "{pdb_align}/split_chain/redo/{pdb_align}_{chain_id}_{uniprot_seq}.fasta",
        "{pdb_align}/uniprot_sequences/{uniprot_seq}.fasta"
     output:
        "{pdb_align}/alignment/redo/{pdb_align}_{chain_id}_{uniprot_seq}.clu"
     run:
        if str(wildcards.uniprot_seq) in str(input[1]):
            shell("cat {wildcards.pdb_align}/split_chain/redo/{wildcards.pdb_align}_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta {wildcards.pdb_align}/uniprot_sequences/{wildcards.uniprot_seq}.fasta > {wildcards.pdb_align}/alignment/redo/input_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta")
            shell("{config[clustlo_bin]}  -i {wildcards.pdb_align}/alignment/redo/input_{wildcards.chain_id}_{wildcards.uniprot_seq}.fasta -o {output} --outfmt=clustal --resno")

         
rule pdb4amber:
    input:
         "{pdb}/{pdb}_original.pdb"
    output:
        "{pdb}/pdb4amber/original/{pdb}_amber_original.pdb"
    shell:
        """
        set +eu
        source /usr/local/amber-20/amber.sh
        set -eu
        pdb4amber -i  "{wildcards.pdb}/{wildcards.pdb}_original.pdb" -o {output} -l {wildcards.pdb}/pdb4amber/original/log_file
        """
rule pdb4amber_redo:
    input:
         "{pdb_redo}/{pdb_redo}_redo.pdb"
    output:
        "{pdb_redo}/pdb4amber/redo/{pdb_redo}_amber_redo.pdb"
    shell:
        """
        set +eu
        source /usr/local/amber-20/amber.sh
        set -eu
        pdb4amber -i  "{wildcards.pdb_redo}/{wildcards.pdb_redo}_redo.pdb" -o {output} -l {wildcards.pdb_redo}/pdb4amber/redo/log_file
        """

