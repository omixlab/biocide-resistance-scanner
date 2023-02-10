import argparse
from argparse import ArgumentParser
import subprocess, os

parser = argparse.ArgumentParser(
                    prog = 'BioBacResist',
                    description = 'A tool to indetify biocide resistant genes in bacteria')
def main ():
    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input FASTA/Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output CSV file')
    argument_parser.add_argument('--abricate-setupdb', help='Runs abricate database setup', action='store_true', default=False)
    arguments = argument_parser.parse_args()

    run_abricate(arguments.input, os.path.join(arguments.output, 'abricate.tab'), arguments.abricate_setupdb)

def run_abricate (input_file, output_file, abricate_setupdb):
    if abricate_setupdb == True:
        subprocess.call(f"abricate --setupdb", shell=True)
    subprocess.call(f"abricate --db plasmidfinder {input_file} > {output_file}", shell=True)


def run_islandpath (input_file, output_file):
    

if __name__ == '__main__':
    main()