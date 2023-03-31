import subprocess
import os

DEVNULL = open(os.devnull)

def run_platon (input_file, output_directory):
    subprocess.call(f'platon --db data/dbs/platon/db/ -o {output_directory} {input_file}', shell=True, stdout=DEVNULL, stderr=DEVNULL)
