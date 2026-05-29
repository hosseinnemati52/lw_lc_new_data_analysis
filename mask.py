# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:32:57 2026

@author: Nemat002
"""

import pandas as pd
import sys
from pathlib import Path
import numpy as np


def excel_to_dict(file_path, sheet_name, result_init={}):
    # Skip the header row
    df = pd.read_excel(file_path, skiprows=1, header=None, sheet_name=sheet_name)

    # result = {}
    result = result_init.copy()

    for _, row in df.iterrows():
        key = row[0]      # name column
        
        if key not in result.keys():
            obj_type = row[1] # m/w/c column
        else:
            print("Error!")
            sys.exit()

        values = [
            int(v)
            for v in row[2:]
            if pd.notna(v)
        ]

        result[key] = {
            "type": obj_type,
            "nokay_frames": values
        }

    return result

def mask_func(read_address, write_address, problems_dict):
    
    def do_for_matrix(matrix):
        new_matrix = matrix.copy()
        n_orgs  = np.shape(matrix)[0]
        n_times = np.shape(matrix)[1]
        if np.max(matrix)>0:
            for org_c in range(n_orgs):
                org_key = org_names[org_c]
                if org_key not in problems_dict.keys():
                    continue
                for t_c in range(n_times):
                    if t_c in problems_dict[org_key]['nokay_frames']:
                        new_matrix[org_c, t_c] = None
        return new_matrix
    
    with open(read_address+"/org_names.txt", "r") as org_names_file:
        org_names = [line.strip() for line in org_names_file]
        
    w_mat = np.loadtxt(read_address+'/w_mat.txt', delimiter=',')
    c_mat = np.loadtxt(read_address+'/c_mat.txt', delimiter=',')
    
    Path(write_address).mkdir(exist_ok=True)
    with open(write_address+"/org_names.txt", "w") as org_names_file:
        for org_name in org_names:
            org_names_file.write(org_name + "\n")
    
    w_mat_new = do_for_matrix(w_mat)
    c_mat_new = do_for_matrix(c_mat)
    
    np.savetxt(write_address+'/w_mat.txt', X=w_mat_new, delimiter=',', fmt='%.0f')
    np.savetxt(write_address+'/c_mat.txt', X=c_mat_new, delimiter=',', fmt='%.0f')

    return

problems_dict = {}
file_path = 'Sebas_new_data_orig/orgs_summary.xlsx'
sheet_name = 'so-so'
problems_dict = excel_to_dict(file_path, sheet_name, problems_dict)

sheet_name = 'problematic'
problems_dict = excel_to_dict(file_path, sheet_name, problems_dict)

read_address = 'Sebas_new_data_orig/pure_w'
write_address = 'Sebas_new_data_filtered/pure_w'
mask_func(read_address, write_address, problems_dict)

read_address = 'Sebas_new_data_orig/pure_c'
write_address = 'Sebas_new_data_filtered/pure_c'
mask_func(read_address, write_address, problems_dict)

read_address = 'Sebas_new_data_orig/mix'
write_address = 'Sebas_new_data_filtered/mix'
mask_func(read_address, write_address, problems_dict)

