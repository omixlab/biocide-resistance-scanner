from brs.lib_ssg_lugia.main import SSG_LUGIA
from Bio import SeqIO
import json
import os

def run_ssg_lugia(input_file, output_directory):

    islands = {}
    for r, record in enumerate(SeqIO.parse(input_file, 'fasta')):
        temp_record_fasta = os.path.join(output_directory, f'genome_{r}.fasta')
        SeqIO.write([record], temp_record_fasta, 'fasta')
        islands[record.id] = SSG_LUGIA(sequence_fasta_file_path=temp_record_fasta, model_name='SSG-LUGIA-F')
    with open(os.path.join(output_directory, 'ssg_lugia.json'), 'w') as writer:
        writer.write(json.dumps(islands))
    return islands
