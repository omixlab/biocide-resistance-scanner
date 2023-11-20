from Bio.Blast import NCBIXML
import pandas as pd
import subprocess
import os

DEVNULL = open(os.devnull)

def run_blast_bacmet(input_fasta_file, output_directory, db='bacmet2_exp', evalue=1e-100, min_identity=0.7):
    
    output_xml_file = os.path.join(output_directory, 'bacmet.xml')
    output_csv_file = os.path.join(output_directory, 'bacmet.csv')

    subprocess.call(f"blastp -query {input_fasta_file} -evalue {evalue} -db data/dbs/bacmet2/{db} -outfmt 5 -out {output_xml_file}", shell=True, stdout=DEVNULL, stderr=DEVNULL)
        
    bacmet_mapping = pd.read_csv('data/dbs/bacmet2/BacMet2_EXP.753.mapping.txt', sep='\t') 

    bacmet_hits = []

    for record in NCBIXML.parse(open(output_xml_file)): 
        for alignment in record.alignments:
            bacmet_hit_data = {'query_id': record.query}
            bacmet_id =  alignment.hit_def.split("|")[0]
            bacmet_hit_data['bacmet_id'] = alignment.hit_def.split("|")[0]
            bacmet_hit_info = bacmet_mapping.query('BacMet_ID == @bacmet_id').iloc[0]
            bacmet_hit_data['gene_name'] = bacmet_hit_info.Gene_name
            bacmet_hit_data['accession'] = bacmet_hit_info.Accession
            bacmet_hit_data['organism'] = bacmet_hit_info.Organism
            bacmet_hit_data['location'] = bacmet_hit_info.Location
            bacmet_hit_data['compound'] = bacmet_hit_info.Compound
            bacmet_hit_data['identity'] = alignment.hsps[0].identities / alignment.hsps[0].align_length
            bacmet_hit_data['e-value'] = alignment.hsps[0].expect

            if bacmet_hit_data['identity'] > 0.7:
                bacmet_hits.append(bacmet_hit_data)
            break

    df_bacmet_hits = pd.DataFrame(bacmet_hits)
    df_bacmet_hits.to_csv(output_csv_file, index=False)
    
    return pd.DataFrame(bacmet_hits)
