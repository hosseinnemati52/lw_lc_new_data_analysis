#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os

lw_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]

# lw_list = [1]
# lc_list = [1]

mother_folder = os.getcwd()
os.chdir(mother_folder)

# pure_betas = np.loadtxt('pure_organoid_exponents.csv', delimiter=',')
# beta_w_0 = pure_betas[0,0]
# # beta_w_0 = 0.0284
# beta_c_0 = pure_betas[0,1]


v_c_avg_w = 1035 # um^3
v_c_avg_c = 919  # um^3
shell_thickness = 10 #um

    
# sorted_names = ["a", "b", "c", "d"]
# unsorted_names = ["c", "a", "d", "b"]

# connected_list = [30, 10, 40, 20]
# matrix = np.array([
#     [3, 3],
#     [1, 1],
#     [4, 4],
#     [2, 2],
# ])

# # order[i] gives the index in unsorted_names of sorted_names[i]
# index = {name: i for i, name in enumerate(unsorted_names)}
# order = [index[name] for name in sorted_names]

# new_names = [unsorted_names[i] for i in order]
# new_connected_list = [connected_list[i] for i in order]
# new_matrix = matrix[order, :]

# print(new_names)
# print(new_connected_list)
# print(new_matrix)

for l_w in lw_list:
    for l_c in lc_list:
        
        # go to daughter folder
        folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        try:
            os.chdir(mother_folder + "\\"+ folder_name) #windows
        except:
            os.chdir(mother_folder + "/"+ folder_name) #linux
        
        time_points = np.loadtxt("time_points.txt", delimiter=',')
        
        # content = []
        with open("org_names.txt", "r") as org_names_file:
            org_names = [line.strip() for line in org_names_file]
        
        sorted_org_names = sorted(set(org_names), key=str.lower)
        
        
        w_mat = np.loadtxt("w_mat.txt", delimiter=',')
        c_mat = np.loadtxt("c_mat.txt", delimiter=',')
        w_a_mat = np.loadtxt("w_a_mat.txt", delimiter=',')
        c_a_mat = np.loadtxt("c_a_mat.txt", delimiter=',')
        w_v_mat = np.loadtxt("w_v_mat.txt", delimiter=',')
        c_v_mat = np.loadtxt("c_v_mat.txt", delimiter=',')
        w_b_mat = np.loadtxt("w_b_mat.txt", delimiter=',')
        c_b_mat = np.loadtxt("c_b_mat.txt", delimiter=',')
        avg_num_w_visib_to_aff_c_mat  = np.loadtxt("avg_num_w_visib_to_aff_c_mat.txt", delimiter=',')
        avg_num_c_visib_to_aff_w_mat  = np.loadtxt("avg_num_c_visib_to_aff_w_mat.txt", delimiter=',')
        
        index = {name: i for i, name in enumerate(org_names)}
        order = [index[name] for name in sorted_org_names]
        
        w_mat = w_mat[order, :]
        c_mat = c_mat[order, :]
        w_a_mat = w_a_mat[order, :]
        c_a_mat = c_a_mat[order, :]
        w_v_mat = w_v_mat[order, :]
        c_v_mat = c_v_mat[order, :]
        w_b_mat = w_b_mat[order, :]
        c_b_mat = c_b_mat[order, :]
        avg_num_w_visib_to_aff_c_mat  = avg_num_w_visib_to_aff_c_mat[order, :]
        avg_num_c_visib_to_aff_w_mat  = avg_num_c_visib_to_aff_w_mat[order, :]

        
        np.savetxt("w_mat.txt", X=w_mat, fmt='%d', delimiter=',')
        np.savetxt("c_mat.txt", X=c_mat, fmt='%d', delimiter=',')
        np.savetxt("w_a_mat.txt", X=w_a_mat, fmt='%d', delimiter=',')
        np.savetxt("c_a_mat.txt", X=c_a_mat, fmt='%d', delimiter=',')
        np.savetxt("w_v_mat.txt", X=w_v_mat, fmt='%d', delimiter=',')
        np.savetxt("c_v_mat.txt", X=c_v_mat, fmt='%d', delimiter=',')
        np.savetxt("w_b_mat.txt", X=w_b_mat, fmt='%d', delimiter=',')
        np.savetxt("c_b_mat.txt", X=c_b_mat, fmt='%d', delimiter=',')
        np.savetxt("avg_num_w_visib_to_aff_c_mat.txt", X=avg_num_w_visib_to_aff_c_mat, fmt='%.3e', delimiter=',')
        np.savetxt("avg_num_c_visib_to_aff_w_mat.txt", X=avg_num_c_visib_to_aff_w_mat, fmt='%.3e', delimiter=',')
        
        with open("org_names.txt", "w") as org_names_file_sorted:
            for org_name in sorted_org_names:
                org_names_file_sorted.write(org_name + "\n")
                
        os.chdir(mother_folder)

        
        
        