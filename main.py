#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:35:05 2024

@author: Ann Wang
"""

import argparse
from Bio import SeqIO, Seq
import pandas as pd
import os
import warnings
from lib.log import logger
import lib.utlis_sequences as utlis_seqs
import config
from lib.utlis_blastp import run_blastp, save_blastp_output
warnings.simplefilter('ignore')


# Structureing useful information into specific columns, and remian others in to 'others_info' columns 
def parsing_info(cds_infos):
    structured_info = dict(GeneID=[], protein=[], ProteinID=[], others_info=[])
    try:
        for infos in cds_infos:
            [structured_info[key].append(None) for key in structured_info.keys()]
            others_info = []
            for info in infos:
                others_info_condition = [info.find("cds") == -1, info.find("location") == -1,
                                          info.find("gbkey=CDS") == -1, info.find("gene=") == -1,
                                          info.find("db_xref") == -1, info.find("protein") == -1]
 
                if info.find("GeneID") != -1:
                    structured_info["GeneID"][-1] = info.split("GeneID:")[-1]

                if info.find("protein=") != -1:
                    structured_info["protein"][-1] = info.split("protein=")[-1]
 
                if info.find("protein_id") != -1:
                    structured_info["ProteinID"][-1] = info.split("protein_id=")[-1]
                
                if all(others_info_condition) :
                    others_info.append(info)

                structured_info["others_info"][-1] = " / ".join(others_info)
    except Exception as e:
        raise str(e)
    return structured_info

"""
Building table with each records' sequences id, complete coding sequences(CDS), 
translated amino acid sequences, and structured informations
"""
def create_info_table(cds_dataset):
    cds_ids = ["_".join(record.id.split("cds_")[-1].split("_")[0:-1]) for record in cds_dataset]
    cds_seqs = [str(record.seq) for record in cds_dataset]
    pp_seqs = [str(Seq.translate(record.seq)) for record in cds_dataset]
    cds_infos = [record.description.replace("]","").split(" [") for record in cds_dataset]
    structured_info = parsing_info(cds_infos)
    if len(cds_ids) == len(structured_info["GeneID"]):
        bigdata = dict(Contig_ID=cds_ids, Contig_seq=cds_seqs, Contig_pp_seq=pp_seqs)
        # merge the sequences and informations
        bigdata.update(structured_info)
    df = pd.DataFrame(bigdata).set_index("Contig_ID")
    return df

# removing partial cds
# including sequences either have error start_codon, error stop_codonm or is not divisible by 3
def partial_cds_detection(df):
    print("Detecting partial sequences...")
    logger.info("Detecting partial sequences...")
    try:
        df["remainder"] = df.Contig_seq.apply(lambda x:len(x)%3)
        df['start_codon'] = df.Contig_seq.apply(lambda x:x[0:3])
        df['end_codon'] = df.Contig_seq.apply(lambda x:x[-3:])
        has_remainder_cds_idx = df[df["remainder"]!=0].index
        if has_remainder_cds_idx.to_list():
            for record in has_remainder_cds_idx:
                logger.info(f"ContigID {record} is not divisible by 3.")
        con = (df['start_codon']!="ATG") | ((df['end_codon'] != "TAA") & (df['end_codon'] != "TAG") & (df['end_codon'] != "TGA"))
        partial_cds_idx = df[con].index
        if partial_cds_idx.to_list():
            for record in partial_cds_idx:
                logger.info(f"ContigID {record} is partial on 5' or 3'.")
        remove_idx = list(set(has_remainder_cds_idx.to_list() + partial_cds_idx.to_list()))
        remove_df = df.loc[remove_idx].drop(["remainder","start_codon","end_codon"], axis=1)
        complete_df = df.drop(remove_idx).drop(["remainder","start_codon","end_codon"], axis=1)
    except Exception as e:
        raise str(e)
    print(f"Done. Input {len(df)} Contigs, {len(remove_df)} had been removed and {len(complete_df)} Contigs remained in dataset.")
    logger.info(f"Done. Input {len(df)} Contigs, {len(remove_df)} had been removed and {len(complete_df)} Contigs remained in dataset.")
    return complete_df, remove_df

# removing dulicates sequences according to coding sequences or translated aa sequences
def remove_duplicates_seqs(df, target_col):
    try:
        if target_col == "Contig_seq":
            print("Removing duplicated sequences...")
            logger.info("Removing duplicated DNA sequences...")
        elif target_col == "Contig_pp_seq":
            print("Removing duplicated peptide sequences...")
            logger.info("Removing duplicated peptide sequences...")
    
        clean_df = pd.DataFrame(columns=df.columns)
        duplicates_df = pd.DataFrame(columns=df.columns.to_list()+["SameSeqContigID"])
        for seq, group in df.groupby(target_col):
            d = group.iloc[0]
            clean_df.loc[d.name] = d
            if len(group) != 1:
                duplicates = group.iloc[1:]
                duplicates["SameSeqContig"] = [d.name]*(len(duplicates))
                duplicates_df = pd.concat([duplicates_df, duplicates])
    except Exception as e:
        raise str(e)            
    print(f"Done. Input {len(df)} Contigs, {len(duplicates_df)} had been removed and {len(clean_df)} Contigs remained in dataset.")
    logger.info(f"Done. Input {len(df)} Contigs, {len(duplicates_df)} had been removed and {len(clean_df)} Contigs remained in dataset.")
    return clean_df, duplicates_df

# saving main clean complete_cds_table, and duplicates sequences into other tables (for future analysis if needed)
def save_contig_info_table(complete_cds_table, duplicates_table, duplicates_pp_table, output_folder):
    try:
        if not os.path.isdir(output_folder):
            os.makedirs(output_folder, exist_ok=True)
    except Exception as e:
        raise str(e)
    try:
        complete_cds_table.to_csv(os.path.join(output_folder,"complete_cds_table.csv"))
        duplicates_table.to_csv(os.path.join(output_folder,"duplicates_table.csv"))
        duplicates_pp_table.to_csv(os.path.join(output_folder,"duplicates_pp_table.csv"))
    except Exception as e:
        raise str(e)

# main function
def cds_datamining(cds_filepath, output_folder):
    filename = cds_filepath.split("/")[-1]
    logger.info(f"Target fasta file: {filename}")
    cds_dataset = list(SeqIO.parse(cds_filepath, format="fasta"))
    df = create_info_table(cds_dataset)
    complete_df, remove_df = partial_cds_detection(df)
    remove_duplicates_df, duplicates_df = remove_duplicates_seqs(complete_df, target_col="Contig_seq")
    remove_isoform_df, isoform_df = remove_duplicates_seqs(remove_duplicates_df, target_col="Contig_pp_seq")
    save_contig_info_table(remove_isoform_df, duplicates_df, isoform_df, output_folder)
    return remove_isoform_df, duplicates_df, isoform_df


if __name__ == "__main__":
    config.reload_config()

    parser = argparse.ArgumentParser()
    parser.add_argument("CDS_filepath", help="Filepath of target CDS, the file format need to be '.fasta' ,'.fna', or '.fsa'.")
    parser.add_argument("Output_folder_name", help="Name a folder for saving all output data.")
    parser.add_argument("Output_fna_filename", help="Name your clean cds fasta file.")
    parser.add_argument("--sql_dir", help="The filepath and name of the sqlite dataset you want to save your blastp outputs.")
    
    
    args = parser.parse_args()
    
    CDS_filepath = args.CDS_filepath
    Output_folder_name = args.Output_folder_name
    Output_fasta_filename = args.Output_fna_filename
    Sql_dir = args.sql_dir
    
    logger.info("=== programming initiated ...===")
    
    if not Sql_dir:
        # step1: cleaning cds
        #clean_df, remove_duplicates_df, isoform_df = cds_datamining(CDS_filepath, Output_folder_name) 
        #utlis_seqs.save_clean_cds_pp_fasta(clean_df, Output_fasta_filename, Output_folder_name)
    
        # step2: executing blastp
        input_fasta = os.path.join(Output_folder_name, f"{Output_fasta_filename}.faa")
        p = run_blastp(input_fasta, Output_folder_name)
        
    
    
    else:
        # step3: save blastp output to sql
        blastp_output = os.path.join(Output_folder_name, f"{Output_fasta_filename}.out")
        save_blastp_output(blastp_output, Sql_dir, Output_fasta_filename)
    
