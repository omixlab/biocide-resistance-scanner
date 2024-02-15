# biocide-resistance-scanner

A tool for detection of biocide resistance genes in horizontally tranfered regions. 

## Istalling dependencies

```
$ git clone https://github.com/omixlab/biocide-resistance-scanner
$ cd biocide-resistance-scanner
$ conda env create --file environment.yml
$ conda activate biocide-resistance-scanner
$ make setup_databases
```
## Using the tool

```
$ conda activate biocide-resistance-scanner
$ brs --help
usage: brs [-h] --input INPUT --output OUTPUT [--evalue EVALUE]
           [--identity IDENTITY]
           [--abricate-database {plasmidfinder,argannot,card,megares,ncbi,vfdb,ecoh,ecoli_vf,resfinder}]
           [--no-abricate] [--no-bacmet] [--no-platon] [--no-ss-lugia]
           [--no-phispy]

Bacterial genome in a FASTA/Genbank file

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Input Genbank file
  --output OUTPUT       Output directory
  --evalue EVALUE       BLAST evalue threshold (default: 1e-100)
  --identity IDENTITY   BLAST identity threshold (default: 0.7)
  --abricate-database {plasmidfinder,argannot,card,megares,ncbi,vfdb,ecoh,ecoli_vf,resfinder}
                        ABRICATE database to be used (default: ncbi)
  --no-abricate         do not run abricate
  --no-bacmet           do not run bacmet
  --no-platon           do not run platon
  --no-ss-lugia         do not run SS-LUGIA
  --no-phispy           do not run Phispy
```

