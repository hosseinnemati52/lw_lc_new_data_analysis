#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 14:15:07 2026

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


mother_folder = os.getcwd()
os.chdir(mother_folder)
# go to daughter folder
folder_name = 'lw_'+str(int(1))+'__lc_'+str(int(1))
try:
    os.chdir(mother_folder + "\\"+ folder_name) #windows
except:
    os.chdir(mother_folder + "/"+ folder_name) #linux
    
with open("org_names.txt", "r") as org_names_file:
    org_names = [line.strip() for line in org_names_file]
    
w_mat = np.loadtxt("w_mat.txt", delimiter=',')
w_b_mat = np.loadtxt("w_b_mat.txt", delimiter=',')
c_mat = np.loadtxt("c_mat.txt", delimiter=',')
c_b_mat = np.loadtxt("c_b_mat.txt", delimiter=',')

os.chdir(mother_folder)


a_w = 100
a_c = 100

for i in range(len(org_names)):

    
    A_w = a_w * w_mat
    A_c = a_c * c_mat
    
    A_w_touch = a_w * w_b_mat
    A_c_touch = a_c * c_b_mat
    
    A_tot = A_w + A_c
    
    r = (A_tot/(4*np.pi))**0.5
    
    th_m = np.arccos( (A_w-A_c) / A_tot)
    
    argument =( (A_w + A_c_touch) \
              - (A_c - A_c_touch) ) \
                / A_tot
                
    th_c = np.arccos( argument )
    
    argument =( (A_w - A_w_touch) \
              - (A_c + A_w_touch) ) \
                / A_tot
                
    th_w = np.arccos( argument )
    
    width_w = ( (th_w - th_m)*r ) / np.sqrt(a_w)
    
    width_c = ( (th_m - th_c)*r ) / np.sqrt(a_c)

x = np.array(np.mean(width_w, axis=1))
y = np.array(np.mean(width_c, axis=1))

bins = np.quantile(x, [0, 1/3, 2/3, 1])
groups = np.digitize(x, bins[1:-1])

list1 = x[groups == 0]
list2 = x[groups == 1]
list3 = x[groups == 2]

plt.scatter(list1, y[groups == 0])
plt.scatter(list2, y[groups == 1])
plt.scatter(list3, y[groups == 2])


low_jg_orgs = []
mid_jg_orgs = []
hii_jg_orgs = []

low_jg_vals = []
mid_jg_vals = []
hii_jg_vals = []

for i in range(len(org_names)):
    
    if groups[i]==0:
        low_jg_orgs.append(org_names[i])
        low_jg_vals.append(x[i])
    elif groups[i]==1:
        mid_jg_orgs.append(org_names[i])
        mid_jg_vals.append(x[i])
    elif groups[i]==2:
        hii_jg_orgs.append(org_names[i])
        hii_jg_vals.append(x[i])
        
        
with open("low_jg_orgs.txt", "w") as f1:
    for org_name in low_jg_orgs:
        f1.write(org_name + "\n")
        
with open("mid_jg_orgs.txt", "w") as f2:
    for org_name in mid_jg_orgs:
        f2.write(org_name + "\n")

with open("hii_jg_orgs.txt", "w") as f3:
    for org_name in hii_jg_orgs:
        f3.write(org_name + "\n")
    
np.savetxt("low_jg_vals.txt", X=low_jg_vals, fmt='%.3e', delimiter=',')
np.savetxt("mid_jg_vals.txt", X=mid_jg_vals, fmt='%.3e', delimiter=',')
np.savetxt("hii_jg_vals.txt", X=hii_jg_vals, fmt='%.3e', delimiter=',')