import subprocess
import os

DEVNULL = open(os.devnull)

def run_abricate (input_file, output_file, database='card'):
    subprocess.call(f"abricate --db {database} {input_file} > {output_file}", shell=True, stdout=DEVNULL, stderr=DEVNULL)
