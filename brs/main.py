from brs.runners import *
from brs.utils import *
from brs.merger import combine_results
import argparse
from argparse import ArgumentParser
import subprocess, os
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import json

DEVNULL = open(os.devnull, 'w')

parser = argparse.ArgumentParser(
                    prog = 'BiocideResistScanner',
                    description = 'A tool to identify biocide resistant genes in bacteria')
def main():

    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output directory')
    argument_parser.add_argument('--evalue', default=1e-100, help='BLAST evalue threshold (default: 1e-100)')
    argument_parser.add_argument('--identity', default=0.7, help='BLAST identity threshold (default: 0.7)')
    argument_parser.add_argument('--abricate-database', default='ncbi', help='ABRICATE database to be used (default: ncbi)', choices=['plasmidfinder', 'argannot', 'card', 'megares', 'ncbi', 'vfdb', 'ecoh', 'ecoli_vf', 'resfinder'])
    argument_parser.add_argument('--no-abricate', default=False, action='store_true', help='do not run abricate')
    argument_parser.add_argument('--no-bacmet', default=False, action='store_true', help='do not run bacmet')
    argument_parser.add_argument('--no-platon', default=False, action='store_true', help='do not run platon')
    argument_parser.add_argument('--no-ss-lugia', default=False, action='store_true', help='do not run SS-LUGIA')
    argument_parser.add_argument('--no-phispy', default=False, action='store_true', help='do not run Phispy')
    arguments = argument_parser.parse_args()

    input_genome_fasta = os.path.join(arguments.output, 'genome.fasta')
    input_proteins_fasta = os.path.join(arguments.output, 'proteins.fasta')

    bacmet_blast_xml = os.path.join(arguments.output, 'bacmet.xml')

    run_converter_gbk_to_genome_fna(arguments.input, input_genome_fasta)
    run_converter_gbk_to_protein_fna(arguments.input, input_proteins_fasta)

    if len(list(SeqIO.parse(input_proteins_fasta, 'fasta'))) == 0:
        print(f'Error: "{arguments.input}" contains no annotated proteins')
        exit(1)

    if not arguments.no_abricate:
        print('Running Abricate to find antibiotic resistance genes ...')
        _  = run_abricate(arguments.input, os.path.join(arguments.output, 'abricate.tab'))
 
    if not arguments.no_bacmet:
        print('Running BLAST against BacMet to find biocide resistance genes ...')
        _ = run_blast_bacmet(input_proteins_fasta, arguments.output, min_identity=arguments.identity, evalue=arguments.evalue)

    if not arguments.no_platon:
        print('Running Platon to find plasmids ...')
        _ = run_platon(input_genome_fasta, arguments.output)

    if not arguments.no_ss_lugia:
        print('Running SSG-LUGIA to find genomic islands ...')
        _ = run_ssg_lugia(arguments.input, arguments.output)

    if not arguments.no_phispy:
        print('Running Phispy to find phages ...')
        _ = run_phispy(arguments.input, os.path.join(arguments.output, 'phispy'))

    print('Merging results ...')

    results = combine_results(arguments.output)
    results.to_csv(os.path.join(arguments.output, 'results.csv'), index=False)



if __name__ == '__main__':
    main()
