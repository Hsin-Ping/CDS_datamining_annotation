#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 17:39:27 2024

@author: wangxinping
"""
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from lib.log import logger

def output_fasta_file(seq_list, filename, filepath, ftype="fna"):
    print(f"Writing fasta file: {filename}...")
    logger.info(f"Writing fasta file: {filename}...")
    if filepath:
        if not os.path.isdir(filepath):
            os.makedirs(filepath, exist_ok=True)
        filename = os.path.join(filepath, filename)
    with open(f"{filename}.{ftype}", "w+") as f:
        SeqIO.write(seq_list, f, "fasta")
    print("Done.")
    logger.info("Done")
    return None


def save_clean_cds_pp_fasta(clean_df, fasta_name, filepath):
    cds_dataset = []
    pp_dataset = []
    try:
        
        for idx in range(len(clean_df)):
            record = clean_df.iloc[idx]
            cds = SeqRecord(record.Contig_seq, id = record.name, description=f"[{record.protein}]")
            pp = SeqRecord(record.Contig_pp_seq, id = record.name, description=f"[{record.protein}]")
            cds_dataset.append(cds)
            pp_dataset.append(pp)
    except Exception as e:
        raise str(e)
    try:
        output_fasta_file(cds_dataset, fasta_name, filepath)
        output_fasta_file(pp_dataset, fasta_name, filepath, ftype="faa")
    except Exception as e:
        raise str(e)