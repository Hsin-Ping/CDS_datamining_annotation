#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 21:24:51 2024

@author: wangxinping
"""
import os
import config
import sqlite3
import pandas as pd
from lib.log import logger
import paramiko
    
def connet_sqlite(sql_name):
    conn = sqlite3.connect(sql_name)
    return conn

def save_blastp_output(sql_name, blastp_output_df, save_sheet_name):
    try:
        conn = connet_sqlite(sql_name)
        blastp_output_df.to_sql(sql_name, conn, if_exists="replace")
        logger.info(f"Successfully save blastp output into sheet '{save_sheet_name}' in sql '{sql_name}'")
    except Exception as e:
        raise str(e)
    finally: 
        conn.close()
        
def connect_docker_server():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(hostname=config.docker_server_ip, port=config.docker_server_port, \
                username=config.docker_server_username, password=config.docker_server_password, \
                key_filename=config.docker_server_key_file)
    _stdin, _stdout,_stder = client.exec_command('cd blast/ ; pwd')
    server_blast_dir = _stdout.read().decode().strip()
    return client, server_blast_dir 

def ftp_input_to_docker_server(client, server_blast_dir, output_folder, input_fasta):
    cmd = f"mkdir {server_blast_dir}/{output_folder}"
    output_fasta = os.path.join(server_blast_dir, output_folder, input_fasta.split("/")[-1])
    client.exec_command(cmd)
    ftp_client=client.open_sftp()
    ftp_client.put(input_fasta, output_fasta)
    ftp_client.close()
    return

def ftp_get_output_from_docker_server(client, server_blast_dir, output_folder, input_fasta):
    ftp_client=client.open_sftp()
    filename = input_fasta.split("/")[-1].replace(".faa", ".out")
    remote_filepath = os.path.join(server_blast_dir, output_folder, filename)
    local_filepath = os.path.join(os.getcwd(), output_folder, filename)
    ftp_client.get(remote_filepath, local_filepath)
    ftp_client.close()
    return 
   
    
    