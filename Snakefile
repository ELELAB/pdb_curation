configfile: "Config.yaml"
import pandas as pd
import requests
from urllib import request as ur
from urllib.error import HTTPError
from pypdb import get_all_info
import pypdb
import os
from Bio.PDB import PDBParser
#bash_command= "source /usr/local/amber-20/amber.sh"
#os.system(bash_command)

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

### to insert in uniprot_sequences rule ###


"""

def pdb_chain_name(input,pdb):
    from Bio.PDB import PDBParser
    x=str(pdb)
    with open(input) as pdb:
        pdb_file=PDBParser().get_structure(x, input)
        l=[]
        for chainn in pdb_file.get_chains():
            i=str(chainn)
            l.append(i[10])
        i=list(uniprot_chain(wildcards.pdb))
        u=dict(zip(i,l))
    return u.values()

"""

#### list to use along with the wildcards in pdb_split_chain rule to assign each list entry to the pdb_redo file #### 
pdb_list = []
ID_entry_redo=[]
for i in pdb_csv['pdb'].str.lower().to_list():
    try:
        ur.urlopen(f"https://pdb-redo.eu/db/{i}/{i}_final.pdb")
    except HTTPError:
        continue
    pdb_list.append(i.upper()+"/"+i.upper()+"_pdbredo.pdb")
    ID_entry_redo.append(i.upper())

### pdb_method table obtaining ###

pdb_csv['method'] = pdb_csv.apply(lambda x: get_all_info(x['pdb'])['exptl'][0]['method'], axis=1)
pdb_csv.to_csv("pdb_method")
list_chain_ID=[]
for i in pdb_csv['pdb'].to_list():
    dictionary=uniprot_chain(i)
    ID=", ".join(f"{k} {v}" for k,v in dictionary.items())
    list_chain_ID.append(ID)    

pdb_csv["chain_uniprot_ID"]=list_chain_ID
pdb_csv.to_csv("pdb_method")

#####

rule all:
    input:
        pdb_list,
        expand("{pdb}/{pdb}_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/split_chain/", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/alignments/", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/pdb4amber/original/{pdb}_amber_original.pdb", pdb=pdb_csv['pdb'].str.upper()),
     #   expand("{pdb_redo}/pdb4amber/redo/{pdb_redo}_pdbredo.pdb", pdb_redo=ID_entry_redo),

rule download_pdb:
    output:
        "{pdb}/{pdb}_original.pdb",
    run:
        with open(f"{output}",'w') as f:
            f.write(pypdb.get_pdb_file(f'{wildcards.pdb}'))

rule download_redo:
    output:
        "{pdb}/{pdb_r}_pdbredo.pdb"
    run:
        ur.urlretrieve(f"https://pdb-redo.eu/db/{str.lower(wildcards.pdb)}/{str.lower(wildcards.pdb)}_final.pdb", output[0])
                     
rule uniprot_sequence:
    input:
        "{pdb}/{pdb}_original.pdb"
    output:
        directory("{pdb}/uniprot_sequences/")
    run:
        shell("mkdir -p {wildcards.pdb}/uniprot_sequences/")
        uniprot_ID=list(uniprot_chain(wildcards.pdb).values())
        pdb_file=PDBParser().get_structure(wildcards.pdb, input[0])
        chain_list=[]
        for chain in pdb_file.get_chains():
            chain_list.append(chain.get_id())
        chain_unip=list(uniprot_chain(wildcards.pdb))
        d=dict(zip(chain_unip,chain_list))
        y=list(d.values())
        for i in range(len(uniprot_ID)):
            ur.urlretrieve(f'https://www.uniprot.org/uniprot/{uniprot_ID[i]}.fasta', f'{wildcards.pdb}/uniprot_sequences/{wildcards.pdb}_{y[i]}.fasta')

rule split_chain_2fasta:
    input:
        "{pdb}/{pdb}_original.pdb",
        expand("{pdb_redo}/{pdb_redo}_pdbredo.pdb", pdb_redo=ID_entry_redo)
    output:
        directory("{pdb}/split_chain")
    shell:
        """
        cd {wildcards.pdb}
        mkdir split_chain
        cd split_chain
        for filename in ../*.pdb; do
            pdb_splitchain $filename
        done
        mkdir original
        if [[ $filename == *"pdbredo"* ]]; then
                 mkdir redo
        fi
        for filename in *.pdb; do
            if [[ $filename == *"pdbredo"* ]]; then
                 sed -i -e '/TER/Q' $filename
                 {config[pdb2fasta_bin]} $filename > redo/$(basename -- "$filename" .pdb).fasta
                 mv $filename redo
            fi
            if [[ $filename == *"original"* ]]; then
                 sed -i -e '/TER/Q' $filename
                 {config[pdb2fasta_bin]} $filename > original/$(basename -- "$filename" .pdb).fasta
                 mv $filename original
            fi
        done
        """ 

rule alignment:
    input:
        expand("{pdb}/uniprot_sequences/", pdb=pdb_csv['pdb'].str.upper()),
        expand("{pdb}/split_chain/", pdb=pdb_csv['pdb'].str.upper()),
    output:
        directory("{pdb}/alignments"),
    shell:
        """
        mkdir -p {wildcards.pdb}/alignments/original
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
                cat {wildcards.pdb}/uniprot_sequences/*_$i.fasta {wildcards.pdb}/split_chain/original/*_$i.fasta > {wildcards.pdb}/alignments/original/input_$i.fasta
            fi
        done 
        for i in {wildcards.pdb}/alignments/original/*.fasta
        do
            n=${{i%.*}}
            n=${{n##*_}}
            {config[clustlo_bin]} -i $i -o {wildcards.pdb}/alignments/original/alignment_original_$n.clu --outfmt=clustal --resno
        done
        if [[ -d {wildcards.pdb}/split_chain/redo ]]; then
             mkdir -p {wildcards.pdb}/alignments/redo
            for i in "${{array[@]}}"
            do
            if  [[ 'grep -q '$i' ${wildcards.pdb}/uniprot_sequences/*.fasta' ]]; then
                cat {wildcards.pdb}/uniprot_sequences/*_$i.fasta {wildcards.pdb}/split_chain/redo/*_$i.fasta > {wildcards.pdb}/alignments/redo/input_$i.fasta
            fi
            done
            for i in {wildcards.pdb}/alignments/redo/*.fasta
            do
                n=${{i%.*}}
                n=${{n##*_}}
                {config[clustlo_bin]} -i $i -o {wildcards.pdb}/alignments/redo/alignment_redo_$n.clu --outfmt=clustal --resno
            done
        fi
        
        """
rule pdb4amber:
    input:
         "{pdb}/{pdb}_original.pdb",
         expand("{pdb_redo}/{pdb_redo}_pdbredo.pdb", pdb_redo=ID_entry_redo)
    output:
        "{pdb}/pdb4amber/original/{pdb}_amber_original.pdb"
    shell:
        """
        set +eu
        source /usr/local/amber-20/amber.sh
        set -eu
        pdb4amber -i  "{wildcards.pdb}/{wildcards.pdb}_original.pdb" -o {output} -l {wildcards.pdb}/pdb4amber/original/log_file
        if [[ -f {wildcards.pdb}/{wildcards.pdb}_pdbredo.pdb ]]; then
            mkdir {wildcards.pdb}/pdb4amber/redo
            pdb4amber -i {wildcards.pdb}/{wildcards.pdb}_pdbredo.pdb -o {wildcards.pdb}/pdb4amber/redo/{wildcards.pdb}_amber_pdbredo.pdb -l {wildcards.pdb}/pdb4amber/redo/log_file
        fi
        """
       
