import argparse
from argparse import ArgumentParser
import subprocess, os
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import os

parser = argparse.ArgumentParser(
                    prog = 'BiocideResistScanner',
                    description = 'A tool to indetify biocide resistant genes in bacteria')
def main ():
    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input FASTA/Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output CSV file')
    argument_parser.add_argument('--abricate-setupdb', help='Runs abricate database setup', action='store_true', default=False)
    arguments = argument_parser.parse_args()

    run_abricate(arguments.input, os.path.join(arguments.output, 'abricate.tab'), arguments.abricate_setupdb)

    input_proteins_fasta = os.path.join(arguments.output, 'proteins.fasta')

    run_converter_gbk_to_fna(arguments.input, input_proteins_fasta)

    blast_xml = os.path.join(arguments.output, 'bacmet2.xml')

    blast_results = run_blast_bacmet(input_proteins_fasta, blast_xml, "bacmet2_exp")
    
    bacmet_mapping = pd.read_csv('data/dbs/bacmet2/BacMet2_EXP.753.mapping.txt', sep='\t') 

    bacmet_hits = []

    for record in blast_results: 
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
    
    pd.DataFrame(bacmet_hits).to_csv(os.path.join(arguments.output, 'bacmet_hits.csv'))




def run_blast_bacmet(input_fasta_file, output_xml_file, db):
    os.system(f"blastp -query {input_fasta_file} -evalue 1e-100 -db data/dbs/bacmet2/{db} -outfmt 5 -out {output_xml_file}")#, shell=True)
    return NCBIXML.parse(open(output_xml_file))

def run_converter_gbk_to_fna (input_gbk_file, output_fna_file):

    output_handle = open(output_fna_file, "w")

    for seq_record in SeqIO.parse(open(input_gbk_file), "genbank") :
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                output_handle.write(">%s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    str(seq_feature.location.extract(seq_record.seq).translate())))
                
    output_handle.close()


def run_abricate (input_file, output_file, abricate_setupdb):
    if abricate_setupdb == True:
        subprocess.call(f"abricate --setupdb", shell=True)
    subprocess.call(f"abricate --db plasmidfinder {input_file} > {output_file}", shell=True)
    

if __name__ == '__main__':
    main()