#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:12:44 2026

@author: hossein
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os
from scipy.stats import pearsonr
import subprocess
import sys

def indices_in_subset(list_1, list_2):
    subset = set(list_2)
    indices = []
    subset_sorted = []
    for i in range(len(list_1)):
        if list_1[i] in subset:
            indices.append(i)
            subset_sorted.append(list_1[i])
            
    return indices, subset_sorted

def filter_by_indices(indices, arr):
    out = np.zeros((len(indices), arr.shape[1]), dtype=arr.dtype)

    for i in range(len(indices)):
        out[i] = arr[indices[i]]

    return out

lw_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]


mother_folder = os.getcwd()
os.chdir(mother_folder)

org_names_filter = []
del org_names_filter
org_names_filter = []
with open("org_names_filter.txt", "r") as f:
    for line in f:
        org_names_filter.append(line.strip())

for l_w in lw_list:
    for l_c in lc_list:
        
        print('lw: '+str(l_w)+', '+'lc: '+str(l_c))
        
        
        # going to the folder for reading
        folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        try:
            os.chdir(mother_folder + "\\"+"orig_data_all"+"\\"+ folder_name) #windows
        except:
            os.chdir(mother_folder + "/"+"orig_data_all"+"/"+ folder_name) #linux
        
        # reading data
        org_names = []
        del org_names
        org_names = []
        with open("org_names.txt", "r") as f:
            for line in f:
                org_names.append(line.strip())
        
        time_points = np.loadtxt("time_points.txt", delimiter=',')
        # n_time = len(time_points)
        
        w_mat = np.loadtxt("w_mat.txt", delimiter=',')
        w_b_mat = np.loadtxt("w_b_mat.txt", delimiter=',')
        w_a_mat = np.loadtxt("w_a_mat.txt", delimiter=',')
        w_v_mat = np.loadtxt("w_v_mat.txt", delimiter=',')
        w_u_mat = np.loadtxt("w_u_mat.txt", delimiter=',')
        
        c_mat = np.loadtxt("c_mat.txt", delimiter=',')
        c_b_mat = np.loadtxt("c_b_mat.txt", delimiter=',')
        c_a_mat = np.loadtxt("c_a_mat.txt", delimiter=',')
        c_v_mat = np.loadtxt("c_v_mat.txt", delimiter=',')
        c_u_mat = np.loadtxt("c_u_mat.txt", delimiter=',')
        
        avg_num_c_visib_to_aff_w_mat = np.loadtxt("avg_num_c_visib_to_aff_w_mat.txt", delimiter=',')
        avg_num_w_visib_to_aff_c_mat = np.loadtxt("avg_num_w_visib_to_aff_c_mat.txt", delimiter=',')
        
        vol_lum = np.loadtxt("vol_lum.txt", delimiter=',')
        
        # fw_beta_a_w = np.loadtxt("fw_beta_a_w.txt", delimiter=',')
        # fc_beta_a_c = np.loadtxt("fc_beta_a_c.txt", delimiter=',')
        
        fw_beta_a_w_integ = np.loadtxt("fw_beta_a_w_integ.txt", delimiter=',')
        fc_beta_a_c_integ = np.loadtxt("fc_beta_a_c_integ.txt", delimiter=',')
        
        Wa_beta_a_w_integ = np.loadtxt("Wa_beta_a_w_integ.txt", delimiter=',')
        Ca_beta_a_c_integ = np.loadtxt("Ca_beta_a_c_integ.txt", delimiter=',')
        # reading data
        
        # back to the mother folder
        os.chdir(mother_folder)
        
        # making the new folders and go into it
        os.makedirs(folder_name, exist_ok=True)
        try:
            os.chdir(mother_folder +"\\"+ folder_name) #windows
        except:
            os.chdir(mother_folder +"/"+ folder_name) #linux
        
        filter_org_indices, org_names_write = indices_in_subset(org_names, org_names_filter)
        with open("org_names.txt", "w") as f:
            for org_name in org_names_write:
                f.write(org_name + "\n")
        
        #filtering
        w_mat   = filter_by_indices(filter_org_indices, w_mat)
        w_b_mat = filter_by_indices(filter_org_indices, w_b_mat)
        w_a_mat = filter_by_indices(filter_org_indices, w_a_mat)
        w_v_mat = filter_by_indices(filter_org_indices, w_v_mat)
        w_u_mat = filter_by_indices(filter_org_indices, w_u_mat)
        c_mat   = filter_by_indices(filter_org_indices, c_mat)
        c_b_mat = filter_by_indices(filter_org_indices, c_b_mat)
        c_a_mat = filter_by_indices(filter_org_indices, c_a_mat)
        c_v_mat = filter_by_indices(filter_org_indices, c_v_mat)
        c_u_mat = filter_by_indices(filter_org_indices, c_u_mat)
        avg_num_c_visib_to_aff_w_mat = filter_by_indices(filter_org_indices, avg_num_c_visib_to_aff_w_mat)
        avg_num_w_visib_to_aff_c_mat = filter_by_indices(filter_org_indices, avg_num_w_visib_to_aff_c_mat)
        vol_lum = filter_by_indices(filter_org_indices, vol_lum)
        fw_beta_a_w_integ = filter_by_indices(filter_org_indices, fw_beta_a_w_integ)
        fc_beta_a_c_integ = filter_by_indices(filter_org_indices, fc_beta_a_c_integ)
        Wa_beta_a_w_integ = filter_by_indices(filter_org_indices, Wa_beta_a_w_integ)
        Ca_beta_a_c_integ = filter_by_indices(filter_org_indices, Ca_beta_a_c_integ)
        #filtering
        
        # saving
        np.savetxt("time_points.txt", time_points, delimiter=',', fmt='%f')
        np.savetxt("w_mat.txt", w_mat, delimiter=',', fmt='%d')
        np.savetxt("w_b_mat.txt", w_b_mat, delimiter=',', fmt='%d')
        np.savetxt("w_a_mat.txt", w_a_mat, delimiter=',', fmt='%d')
        np.savetxt("w_v_mat.txt", w_v_mat, delimiter=',', fmt='%d')
        np.savetxt("w_u_mat.txt", w_u_mat, delimiter=',', fmt='%d')
        np.savetxt("c_mat.txt", c_mat, delimiter=',', fmt='%d')
        np.savetxt("c_b_mat.txt", c_b_mat, delimiter=',', fmt='%d')
        np.savetxt("c_a_mat.txt", c_a_mat, delimiter=',', fmt='%d')
        np.savetxt("c_v_mat.txt", c_v_mat, delimiter=',', fmt='%d')
        np.savetxt("c_u_mat.txt", c_u_mat, delimiter=',', fmt='%d')
        np.savetxt("avg_num_c_visib_to_aff_w_mat.txt", avg_num_c_visib_to_aff_w_mat, delimiter=',', fmt='%.4e')
        np.savetxt("avg_num_w_visib_to_aff_c_mat.txt", avg_num_w_visib_to_aff_c_mat, delimiter=',', fmt='%.4e')
        np.savetxt("vol_lum.txt", vol_lum, delimiter=',', fmt='%.4e')
        np.savetxt("fw_beta_a_w_integ.txt", fw_beta_a_w_integ, delimiter=',', fmt='%.4e')
        np.savetxt("fc_beta_a_c_integ.txt", fc_beta_a_c_integ, delimiter=',', fmt='%.4e')
        np.savetxt("Wa_beta_a_w_integ.txt", Wa_beta_a_w_integ, delimiter=',', fmt='%.4e')
        np.savetxt("Ca_beta_a_c_integ.txt", Ca_beta_a_c_integ, delimiter=',', fmt='%.4e')
        # saving

        # back to the mother folder
        os.chdir(mother_folder)
    
subprocess.run([sys.executable, "sorter.py"])