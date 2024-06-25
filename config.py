from configparser import ConfigParser
from lib.log import logger
import os

config_file = './config.ini'
#share_config_file = './share/config.ini'

config = ConfigParser()
config.read(config_file)

#share_config = init_share_config()

version = ''

blastp_using_docker = None
blastp_db_folder = ""
blastp_db_name = ""
blastp_max_target_seqs = None
blastp_num_threads = None
blastp_evalue = None
blastp_background = None
blastp_taxids = ""
blastp_remote = None

docker_server_ip = ""
docker_server_username = ""
docker_server_password = ""
docker_server_port = ""
docker_server_key_file = ""

def get_version():
    return version

def check_config_section():
    if not config.has_section('common'):
        config.add_section('common')

    if not config.has_section('path'):
        config.add_section('path')

    if not config.has_section('param'):
        config.add_section('param')
        
    if not config.has_section('blastp'): 
        config.add_section('blastp')

    config.write(open(config_file, 'w'))

def get_config():
    global version, blastp_using_docker, blastp_db_folder, blastp_db_name, blastp_max_target_seqs, blastp_num_threads, blastp_evalue, blastp_background, blastp_taxids, \
        blastp_remote, docker_server_ip, docker_server_username, docker_server_password, docker_server_port, docker_server_key_file
        

    try:
        version = config.get('common', 'version')
    except Exception as e:
        logger.warning(e)
        version = ''
        
    try:
        blastp_using_docker = config.getboolean("blastp", "using_docker")
    except Exception as e:
        logger.warning(e)
        blastp_using_docker = ""
        
    try:
        blastp_db_folder = config.get("blastp", "db_folder")
    except Exception as e:
        logger.warning(e)
        blastp_db_folder = ""
        
    try:
        blastp_db_name = config.get("blastp", "db_name")
    except Exception as e:
        logger.warning(e)
        blastp_db_name = ""    

    try:
        blastp_max_target_seqs = config.get('blastp', 'max_target_seqs')
    except Exception as e:
        logger.warning(e)
        blastp_max_target_seqs = None

    try:
        blastp_num_threads = config.get('blastp', 'num_threads')
    except Exception as e:
        logger.warning(e)
        blastp_num_threads = None
        
    try:
        blastp_evalue = config.get('blastp', 'evalue')
    except Exception as e:
        logger.warning(e)
        blastp_evalue = None
        
    try:
        blastp_background = config.getboolean('blastp', 'background')
    except Exception as e:
        logger.warning(e)
        blastp_background = None
        
    try:
        blastp_taxids = config.get('blastp', 'taxids')
    except Exception as e:
        logger.warning(e)
        blastp_taxids = ""
        
    try:
        blastp_remote = config.get('blastp', 'remote')
    except Exception as e:
        logger.warning(e)
        blastp_remote = ""
        
    try:
        docker_server_ip = config.get('docker_server', 'ip')
    except Exception as e:
        logger.warning(e)
        docker_server_ip = ""
        
    try:
        docker_server_username = config.get('docker_server', 'username')
    except Exception as e:
        logger.warning(e)
        docker_server_username = ""
        
    try:
        docker_server_password = config.get('docker_server', 'password')
    except Exception as e:
        logger.warning(e)
        docker_server_password = ""
        
    try:
        docker_server_port = config.get('docker_server', 'port')
    except Exception as e:
        logger.warning(e)
        docker_server_port = ""

    try:
        docker_server_key_file = config.get('docker_server', 'key_file')
    except Exception as e:
        logger.warning(e)
        docker_server_key_file = ""
        
def reload_config():
    check_config_section()
    get_config()
    #check_config_valid()

if __name__ == '__main__':
    reload_config()
