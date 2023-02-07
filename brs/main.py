import argparse
from argparse import ArgumentParser


parser = argparse.ArgumentParser(
                    prog = 'BioBacResist',
                    description = 'A tool to indetify biocide resistant genes in bacteria',)
def main ():
    argument_parser = ArgumentParser(description='Bacterial genome in a FASTA/Genbank file')
    argument_parser.add_argument('--input', required=True, help='Input FASTA/Genbank file')
    argument_parser.add_argument('--output', required=True, help='Output CSV file')

    arguments = argument_parser.parse_args()

if __name__ == '__main__':
    main()