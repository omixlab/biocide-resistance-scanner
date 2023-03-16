import argparse
from argparse import ArgumentParser
from ssg_lugia.main import SSG_LUGIA
import subprocess, os
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import json

DEVNULL = open(os.devnull, 'w')

parser = argparse.ArgumentParser(
                    prog = 'BiocideResistScanner',
                    description = 'A tool to indetify biocide resistant genes in bacteria')
def main ():

    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input FASTA/Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output CSV file')
    arguments = argument_parser.parse_args()

    input_genome_fasta = os.path.join(arguments.output, 'genome.fasta')
    input_proteins_fasta = os.path.join(arguments.output, 'proteins.fasta')

    bacmet_blast_xml = os.path.join(arguments.output, 'bacmet.xml')

    run_converter_gbk_to_genome_fna(arguments.input, input_genome_fasta)
    run_converter_gbk_to_protein_fna(arguments.input, input_proteins_fasta)

    print('Running Abricate to find antibiotic resistance genes ...')
    abricate_results  = run_abricate(arguments.input, os.path.join(arguments.output, 'abricate.tab'))

    print('Running BLAST against BacMet to find biocide resistance genes ...')
    bacmet_results    = run_blast_bacmet(input_genome_fasta, arguments.output)

    print('Running Platon to find plasmids ...')
    platon_results    = run_platon(input_genome_fasta, arguments.output)
    
    print('Running SSG-LUGIA to find genomic islands ...')
    ssg_lugia_results = run_ssg_lugia(input_genome_fasta, arguments.output)

    print(ssg_lugia_results)

def run_blast_bacmet(input_fasta_file, output_directory, db='bacmet2_exp', evalue=1e-100):
    
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
            bacmet_hit_info['identity'] = alignment.hsps[0].identities / alignment.hsps[0].align_length
            bacmet_hit_info['e-value'] = alignment.hsps[0].expect
            bacmet_hits.append(bacmet_hit_data)
            break

    df_bacmet_hits = pd.DataFrame(bacmet_hits)
    df_bacmet_hits.to_csv(output_csv_file, index=False)
    
    return pd.DataFrame(bacmet_hits)

def run_converter_gbk_to_genome_fna(input_gbk_file, output_fna_file):

    records = SeqIO.parse(input_gbk_file, 'genbank')
    SeqIO.write(records, output_fna_file, 'fasta')    

def run_converter_gbk_to_protein_fna (input_gbk_file, output_fna_file):

    output_handle = open(output_fna_file, "w")

    for seq_record in SeqIO.parse(open(input_gbk_file), "genbank"):
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                output_handle.write(">%s %s:%s:%s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_record.id,
                    str(seq_feature.location.start+1),
                    str(seq_feature.location.end),
                    str(seq_feature.location.extract(seq_record.seq).translate())))
                
    output_handle.close()

def run_abricate (input_file, output_file, database='card'):
    subprocess.call(f"abricate --db {database} {input_file} > {output_file}", shell=True, stdout=DEVNULL, stderr=DEVNULL)

def run_platon (input_file, output_directory):

    subprocess.call(f'platon --db data/dbs/platon/db/ -o {output_directory} {input_file}', shell=True, stdout=DEVNULL, stderr=DEVNULL)

def run_ssg_lugia(input_file, output_directory):

    islands = {}
    for r, record in enumerate(SeqIO.parse(input_file, 'fasta')):
        temp_record_fasta = os.path.join(output_directory, f'genome_{r}.fasta')
        SeqIO.write([record], temp_record_fasta, 'fasta')
        islands[record.id] = SSG_LUGIA(sequence_fasta_file_path=temp_record_fasta, model_name='SSG-LUGIA-F')
    with open(os.path.join(output_directory, 'ssg_lugia.json'), 'w') as writer:
        writer.write(json.dumps(islands))
    return islands

if __name__ == '__main__':
    main()