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
def main ():

    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output CSV file')
    arguments = argument_parser.parse_args()

    input_genome_fasta = os.path.join(arguments.output, 'genome.fasta')
    input_proteins_fasta = os.path.join(arguments.output, 'proteins.fasta')

    bacmet_blast_xml = os.path.join(arguments.output, 'bacmet.xml')

    run_converter_gbk_to_genome_fna(arguments.input, input_genome_fasta)
    run_converter_gbk_to_protein_fna(arguments.input, input_proteins_fasta)

    print('Running Abricate to find antibiotic resistance genes ...')
    _  = run_abricate(arguments.input, os.path.join(arguments.output, 'abricate.tab'))

    print('Running BLAST against BacMet to find biocide resistance genes ...')
    _ = run_blast_bacmet(input_proteins_fasta, arguments.output)

    print('Running Platon to find plasmids ...')
    _ = run_platon(input_genome_fasta, arguments.output)

    print('Running SSG-LUGIA to find genomic islands ...')
    _ = run_ssg_lugia(arguments.input, arguments.output)

    print('Running Phispy to find phages ...')
    _ = run_phispy(arguments.input, os.path.join(arguments.output, 'phispy'))

    print('Merging results ...')
    results = combine_results(arguments.output)
    results.to_csv(os.path.join(arguments.output, 'results.csv'), index=False)





if __name__ == '__main__':
    main()