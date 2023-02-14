import argparse
from argparse import ArgumentParser
import subprocess, os
from Bio import SeqIO
from Bio.Blast import NCBIXML


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

    run_blast_bacmet(input_proteins_fasta, blast_xml)


def run_blast_bacmet(input_fasta_file, output_xml_file):
    subprocess.call(f"blastp -query {input_fasta_file} -db data/dbs/bacmet2/bacmet2_pred -outfmt 5 -out {output_xml_file}", shell=True)
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