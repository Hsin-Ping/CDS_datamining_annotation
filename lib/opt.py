#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 21:24:51 2024

@author: wangxinping
"""

import sqlite3
import pandas as pd
from lib.log import logger

#class sql_dataset:
    
#    def __init__(self, sqlite_name, sheet_name):
#        self.sql = connet_sqlite(sqlite_name)
        
        
    
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