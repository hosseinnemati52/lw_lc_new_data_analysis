# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:48:06 2026

@author: Nemat002
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os

lw_list = [1,2,3,4,5,6,7, 20]
lc_list = [1,2,3,4,5,6,7, 20]

mother_folder = os.getcwd()
os.chdir(mother_folder)

beta_w_0 = 0.0284
beta_c_0 = 0.0398

lw_mat = np.zeros((len(lw_list), len(lc_list)))
lc_mat = np.zeros((len(lw_list), len(lc_list)))

cor_W_Cv   = np.zeros((len(lw_list), len(lc_list)))
cor_W_CvWa = np.zeros((len(lw_list), len(lc_list)))
cor_W_1Wa  = np.zeros((len(lw_list), len(lc_list)))
cor_W_Cv0 = np.zeros((len(lw_list), len(lc_list)))
cor_W_Cv0Wa0 =  np.zeros((len(lw_list), len(lc_list)))
cor_W_nd   = np.zeros((len(lw_list), len(lc_list)))

cor_C_Wv   = np.zeros((len(lw_list), len(lc_list)))
cor_C_WvCa = np.zeros((len(lw_list), len(lc_list)))
cor_C_1Ca  = np.zeros((len(lw_list), len(lc_list)))
cor_C_Wv0 = np.zeros((len(lw_list), len(lc_list)))
cor_C_Wv0Ca0 =  np.zeros((len(lw_list), len(lc_list)))
cor_C_nd   = np.zeros((len(lw_list), len(lc_list)))


for l_w_idx in range(len(lw_list)):
    for l_c_idx in range(len(lc_list)):
        
        l_w = lw_list[l_w_idx]
        l_c = lc_list[l_c_idx]
        
        folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        
        corr_list = np.loadtxt(folder_name+'/corr_list.txt', delimiter=',')
        
        # lw_lc = np.loadtxt(folder_name+'/lw_lc.txt', delimiter=',')
        
        cor_W_Cv[l_w_idx, l_c_idx]   = corr_list[0]
        cor_W_CvWa[l_w_idx, l_c_idx] = corr_list[1]
        cor_W_1Wa[l_w_idx, l_c_idx]  = corr_list[2]
        cor_W_Cv0[l_w_idx, l_c_idx]  = corr_list[3]
        cor_W_Cv0Wa0[l_w_idx, l_c_idx] = corr_list[4]
        cor_W_nd[l_w_idx, l_c_idx]   = corr_list[5]
        
        cor_C_Wv[l_w_idx, l_c_idx] = corr_list[7]
        cor_C_WvCa[l_w_idx, l_c_idx] = corr_list[8]
        cor_C_1Ca[l_w_idx, l_c_idx] = corr_list[9]
        cor_C_Wv0[l_w_idx, l_c_idx] = corr_list[10]
        cor_C_Wv0Ca0[l_w_idx, l_c_idx] = corr_list[11]
        cor_C_nd[l_w_idx, l_c_idx]   = corr_list[12]
        
plt.figure()
plt.plot(lw_list, cor_W_Cv[:,0], label='cor_W_Cv')
plt.plot(lw_list, cor_W_CvWa[:,0], label='cor_W_CvWa')
plt.plot(lw_list, cor_W_1Wa[:,0], label='cor_W_1Wa')
plt.plot(lw_list, cor_W_Cv0[:,0], label='cor_W_Cv0')
plt.plot(lw_list, cor_W_Cv0Wa0[:,0], label='cor_W_Cv0Wa0')
plt.plot(lw_list, cor_W_nd[:,0], label='cor_W_nd')
plt.xlabel('l_w')
plt.ylabel('corr')
plt.legend()
plt.savefig('W_corr.png', dpi=300)


plt.figure()
plt.plot(lc_list, cor_C_Wv[0,:], label='cor_C_Wv')
plt.plot(lc_list, cor_C_WvCa[0,:], label='cor_C_WvCa')
plt.plot(lc_list, cor_C_1Ca[0,:], label='cor_C_1Ca')
plt.plot(lc_list, cor_C_Wv0[0,:], label='cor_C_Wv0')
plt.plot(lc_list, cor_C_Wv0Ca0[0,:], label='cor_C_Wv0Ca0')
plt.plot(lw_list, cor_C_nd[0,:], label='cor_C_nd')
plt.xlabel('l_c')
plt.ylabel('corr')
plt.legend()
plt.savefig('C_corr.png', dpi=300)
