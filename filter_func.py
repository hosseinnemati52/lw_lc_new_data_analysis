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

def indices_in_subset(list_1, list_2):
    subset = set(list_2)
    return [i for i, x in enumerate(list_1) if x in subset]

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
        
        
        
        # back to the mother folder
        os.chdir(mother_folder)
        
        
        
        
        