configfile: "Config.yaml"
import pandas as pd
import requests
from urllib import request as ur
from urllib.error import HTTPError
from pypdb import get_all_info
import pypdb
import os

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
            asym_id
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
            for instance in entity['polymer_entity_instances']:
                final_correspondence[instance['rcsb_polymer_entity_instance_container_identifiers']['asym_id']] = entity['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_accession']
    for i in final_correspondence:
        return final_correspondence
###

pdb_list = []
for i in pdb_csv['pdb'].str.lower().to_list():
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i.upper()+"_pdbredo.pdb")

### pdb_method table obtaining ###

pdb_csv['method'] = pdb_csv.apply(lambda x: get_all_info(x['pdb'])['exptl'][0]['method'], axis=1)
pdb_csv.to_csv("pdb_method")
dictionary:{}
list_chain_ID=[]
for i in pdb_csv['pdb'].to_list():
    dictionary=uniprot_chain(i)
    a=str(dictionary)
    b=a.replace('{', '')
    c=b.replace('}', '')
    string=c.replace("'", '')
    list_chain_ID.append(string)

pdb_csv["chain_uniprot_ID"]=list_chain_ID
pdb_csv.to_csv("pdb_method")

#####

rule all:
    input:
        pdb_list,
        expand("{pdb}/{pdb}_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/alignments/", pdb=pdb_csv['pdb'].str.upper())

rule download_pdb:
    output:
        "{pdb}/{pdb}_original.pdb"
    run:
        with open(f"{output}",'w') as f:
            f.write(pypdb.get_pdb_file(f'{wildcards.pdb}'))
                      

rule download_redo:
    output:
        "{pdb}/{pdb_r}_pdbredo.pdb"
    run:
        ur.urlretrieve(f"https://pdb-redo.eu/db/{str.lower(wildcards.pdb)}/{str.lower(wildcards.pdb)}_final.pdb", output[0])

rule split_chain_2fasta:
    input:
        "{pdb}/{pdb}_original.pdb"
    output:
        directory("{pdb}/split_chain/")
    shell:
        """
        cd {wildcards.pdb}
        mkdir split_chain
        cd split_chain
        pdb_splitchain ../{wildcards.pdb}_original.pdb
        for filename in *.pdb; do
            {config[pdb2fasta_bin]} $filename > $(basename -- "$filename" .pdb).fasta
        done
        """
rule uniprot_sequence:
    output:
        directory("{pdb}/uniprot_sequences/")
    run:
        shell("mkdir -p {wildcards.pdb}/uniprot_sequences/") 
        
        x=list(uniprot_chain(wildcards.pdb).values())
        y=list(uniprot_chain(wildcards.pdb))
        for i in range(len(x)):
            ur.urlretrieve(f'https://www.uniprot.org/uniprot/{x[i]}.fasta', f'{wildcards.pdb}/uniprot_sequences/{wildcards.pdb}_{y[i]}.fasta')

rule alignment:
    input:
        expand("{pdb}/uniprot_sequences/", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/split_chain/", pdb=pdb_csv['pdb'].str.upper()),
    output:
        directory("{pdb}/alignments/")
    shell:
        """
        mkdir -p {wildcards.pdb}/alignments
        declare -a array=()
        for i in {wildcards.pdb}/uniprot_sequences/*.fasta
        do
            n=${{i%.*}}
            n=${{n##*_}}
            array+=("$n")
        done
        for i in "${{array[@]}}"
        do
            if  [[ 'grep -q '$i' ${wildcards.pdb}/uniprot_sequences/*.fasta' ]]; then
                cat {wildcards.pdb}/uniprot_sequences/*_$i.fasta {wildcards.pdb}/split_chain/*_$i.fasta > {wildcards.pdb}/alignments/input_$i.fasta
            fi
        done 
        for i in {wildcards.pdb}/alignments/*.fasta
        do
            n=${{i%.*}}
            n=${{n##*_}}
            {config[clustlo_bin]} -i $i -o {wildcards.pdb}/alignments/alignment_$n.clu --outfmt=clustal --resno
        done
        """
