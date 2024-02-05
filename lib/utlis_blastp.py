#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 18:56:06 2024

@author: wangxinping
"""
import os
import time
import pandas as pd
import config
from lib.log import logger
from Bio import SeqIO
import sqlite3

def run_blastp(input_fasta, output_folder):
    try:
        if not config.blastp_using_docker:
            print("using local software")
            logger.info("run_blastp on local")
            cmd = f"blastp -db {config.blastp_db_folder}/{config.blastp_db_name} -query {input_fasta} " \
                  f"-num_threads {config.blastp_num_threads} -evalue {config.blastp_evalue} " \
                  f"-max_target_seqs {config.blastp_max_target_seqs} " \
                  f"-outfmt '7 qacc sacc stitle pident alignment length qstart qend sstart send evalue bitscore staxids' "\
                  f"-out {output_folder}/{input_fasta.split('/')[-1].split('_pp.fna')[0]}.out"
   
        else:
            print("using docker")
            logger.info("run blastp with docker")
            in_docker_db_folder = config.blastp_db_folder.split("/")[-1]
            in_docker_output_folder = output_folder.split("/")[-1]
            query_name = input_fasta.split("/")[-1]
            cmd = f"sudo docker run --rm -v {config.blastp_db_folder}:/blast/{in_docker_db_folder}:ro" \
                  f" -v {os.getcwd()}/{output_folder}:/blast/{in_docker_output_folder}:rw" \
                  f" ncbi/blast:2.14.0 blastp" \
                  f" -db /blast/{in_docker_db_folder}/{config.blastp_db_name} -query /blast/{in_docker_output_folder}/{query_name}" \
                  f" -evalue {config.blastp_evalue} -max_target_seqs {config.blastp_max_target_seqs} -outfmt '7 qacc sacc stitle pident alignment length qstart qend sstart send evalue bitscore staxids'" \
                  f" -num_threads {config.blastp_num_threads} -out /blast/{in_docker_output_folder}/{input_fasta.split('/')[-1].replace('_pp.fna','.out')}"
        if not config.blastp_taxids == "":
            cmd  = cmd + f" -taxids {config.blastp_taxids}"
        if config.blastp_background:
            cmd = cmd + " &"
    except Exception as e:
        raise str(e)
    try:
        start = time.time()
        print(cmd)
        os.system(cmd)
    except Exception as e:
        logger.error(f"Errors occureed when running blastp...{str(e)}")
    finally:
        t = time.time()-start
        logger.info(f"Done. total time {t}.")
    return None

def get_blastp_output_fields(blastp_output_filepath):
    with open(blastp_output_filepath) as f:
        d = f.readlines()
    fields = list(filter(lambda x:x.find("Fields")!=-1, d))[0]
    fields_name = fields.split(": ")[-1].replace("\n","").split(", ")
    fields_name = list(map(lambda x:x.replace(" ","_").replace(".","").replace("%_",""), fields_name))
    return fields_name
    
def blastp_output_parse(blastp_output_filepath):
    with open(blastp_output_filepath) as f:
        d = f.readlines()
    records = list(filter(lambda x:x.find("#")==-1, d))
    records = list(map(lambda x:x.replace("\n","").split("\t"), records))
    return records


def get_all_blastp_output(blastp_outputs_filepaths):
    records_all = []
    for idx, output in enumerate(blastp_outputs_filepaths):
        if idx==0:
            fields_name = get_blastp_output_fields(output)
            print(fields_name)
        records = blastp_output_parse(output)
        records_all += records
    df = pd.DataFrame(records_all, columns=fields_name)
    return df

def split_fast_file(fasta_file, num_process):
    seqs_list = list(SeqIO.parse(fasta_file, format="fasta"))
    num_seqs_per_process = int(len(seqs_list)/num_process)
    for idx, i in  enumerate(range(0, len(seqs_list),num_seqs_per_process)):
        if idx == num_process-1:
            subset = seqs_list[i:]
        else:
            subset = seqs_list[i:i+num_seqs_per_process]
        print(len(subset))
        with open(f"{fasta_file.split('.fna')[0]}.{idx}.fna", "w+") as f:
            SeqIO.write(subset, f, "fasta")
    return None

def save_blastp_output_to_sql(sql_name, sheetname, blastp_output_df):
    try:
        con = sqlite3.connect(f"{sql_name}.sqlite")
        blastp_output_df.to_sql(sheetname, con, if_exists="replace")
    except:
        print("error")
    finally:
        con.close()
        
def save_blastp_output(blastp_output, sql_name, sheetname):
    try:
        fields_name = get_blastp_output_fields(blastp_output)
        records = blastp_output_parse(blastp_output)
        df = pd.DataFrame(records, columns=fields_name)
        save_blastp_output_to_sql(sql_name, sheetname, df)
    except Exception as e:
        raise str(e)
    return df
