#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:46:32 2024

@author: wangxinping
"""
import os
import paramiko
ssh = paramiko.SSHClient()
ip = "152.70.83.101"
username = "ubuntu"
password = "hsinping"
port = "22"
key_file = "/Users/wangxinping/.ssh/oracle_vm"
#这行代码的作用是允许连接不在know_hosts文件中的主机。
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(hostname=ip, port=port, username=username, password=password, key_filename=key_file)
_stdin, _stdout,_stder = ssh.exec_command('pwd')
server_root_dir = _stdout.read().decode().strip()

input_folder = "Ls_GCF_002870075"
input_filename = "Ls_GCF_002870075"
input_filepath = os.path.join(os.getcwd(), input_folder, f"{input_filename}.faa")
print(input_filepath)
output_filepath = os.path.join(server_root_dir, f"{input_filename}.faa")
print(output_filepath)

ftp_client=ssh.open_sftp()
ftp_client.put(input_filepath, output_filepath)
ftp_client.close()