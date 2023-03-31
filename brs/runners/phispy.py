import subprocess
import os

def run_phispy(input_genbank_file, output_directory):
    os.system(f'mkdir -p {output_directory}')
    subprocess.call(f'phispy  {input_genbank_file} -o {output_directory}', shell=True)