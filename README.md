# CDS_datamining_annotation
This project aim to clean and parse coding sequences from genomic (cds) with any species genome that available in [NCBI genome dataset](https://ftp.ncbi.nlm.nih.gov/genomes/), then using clean, treanslated amino acid sequences to do gene annotation with [pre-formatted blast database](https://ftp.ncbi.nlm.nih.gov/blast/db/) (e.g. landmark database is recommended) by [blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins).<br/>

The procedure includes three steps:
- coding sequences cleaning, structuring, and  translation
- gene annotation with blastp by installing [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or using [ncbi/blast](https://hub.docker.com/r/ncbi/blast/tags) docker image by executing docker command on ubuntu system.
- save gene annotation output in sqlite

```
python3 main.py -h
usage: main.py [-h] [--sql_dir SQL_DIR] CDS_filepath Output_folder_name Output_fna_filename

positional arguments:
  CDS_filepath         Filepath of target CDS, the file format need to be '.fasta' ,'.fna', or '.fsa'.
  Output_folder_name   Name a folder for saving all output data.
  Output_fna_filename  Name your clean cds fasta file.

options:
  -h, --help           show this help message and exit
  --sql_dir SQL_DIR    The filepath and name of the sqlite dataset you want to save your blastp outputs.
```

## Step1: coding sequences cleaning, structuring, and  translation
- preparation: your own genome dataset or download cds from target species genome in [NCBI genome dataset](https://ftp.ncbi.nlm.nih.gov/genomes/)
- workflow:

  <img src="images/step1.png" width=800, height=400></img>


## Step2: gene annotation with blast (Basic Local Alignment Search Tool)
- preparation
  - clean amino acid fasta file (.faa file in first step)
  - installing [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or using [ncbi/blast](https://hub.docker.com/r/ncbi/blast/tags) docker image
  - download the preformatted database provided by NCBI or generate customized database by blastdbcmd command in [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- workflow
## Step3: save gene annotation output in sqlite
