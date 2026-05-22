#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:26:11 2026

@author: hossein
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os

lw_list = [1,2,3,4, 40]
lc_list = [1,2,3,4, 8, 40]

# lw_list = []
# lc_list = [1]

# lw_list = [1]
# lc_list = [1]

mother_folder = os.getcwd()
os.chdir(mother_folder)

pure_betas = np.loadtxt('pure_organoid_exponents.csv', delimiter=',')
beta_w_0 = pure_betas[0,0]
# beta_w_0 = 0.0284
beta_c_0 = pure_betas[0,1]


v_c_avg_w = 1035 # um^3
v_c_avg_c = 919  # um^3
shell_thickness = 10 #um


time_points = np.loadtxt("lw_1__lc_1/time_points.txt", delimiter=',')
n_times = len(time_points)
n_fit_points = 2000

f_beta_w_integ_avg_mat = np.zeros((len(lw_list), n_times))
f_beta_w_integ_err_mat = np.zeros((len(lw_list), n_times))
N_beta_w_integ_avg_mat = np.zeros((len(lw_list), n_times))
N_beta_w_integ_err_mat = np.zeros((len(lw_list), n_times))

w_a_avg_all_lw = np.zeros((len(lw_list), n_times))
w_a_err_all_lw = np.zeros((len(lw_list), n_times))
fw_avg_all_lw = np.zeros((len(lw_list), n_times))
fw_err_all_lw = np.zeros((len(lw_list), n_times))

f_beta_c_integ_avg_mat = np.zeros((len(lc_list), n_times))
f_beta_c_integ_err_mat = np.zeros((len(lc_list), n_times))
N_beta_c_integ_avg_mat = np.zeros((len(lc_list), n_times))
N_beta_c_integ_err_mat = np.zeros((len(lc_list), n_times))

c_a_avg_all_lc = np.zeros((len(lc_list), n_times))
c_a_err_all_lc = np.zeros((len(lc_list), n_times))
fc_avg_all_lc = np.zeros((len(lc_list), n_times))
fc_err_all_lc = np.zeros((len(lc_list), n_times))

F_w_avg = np.zeros((len(lw_list), n_fit_points))
G_w_avg = np.zeros((len(lw_list), n_fit_points))
F_w_err = np.zeros((len(lw_list), n_fit_points))
G_w_err = np.zeros((len(lw_list), n_fit_points))

F_c_avg = np.zeros((len(lc_list), n_fit_points))
G_c_avg = np.zeros((len(lc_list), n_fit_points))
F_c_err = np.zeros((len(lc_list), n_fit_points))
G_c_err = np.zeros((len(lc_list), n_fit_points))

beta_w_a_avg_by_f = np.zeros((len(lw_list), n_fit_points))
beta_w_a_err_by_f = np.zeros((len(lw_list), n_fit_points))
beta_w_a_avg_by_N = np.zeros((len(lw_list), n_fit_points))
beta_w_a_err_by_N = np.zeros((len(lw_list), n_fit_points))

beta_c_a_avg_by_f = np.zeros((len(lc_list), n_fit_points))
beta_c_a_err_by_f = np.zeros((len(lc_list), n_fit_points))
beta_c_a_avg_by_N = np.zeros((len(lc_list), n_fit_points))
beta_c_a_err_by_N = np.zeros((len(lc_list), n_fit_points))

fit_plot_x = np.linspace(np.min(time_points), np.max(time_points), n_fit_points)

# lw_checked = []
# lc_checked = []

outlier_names_list = ['E-139_S-18_O-1']


# W data reading #
for l_w_idx in range(len(lw_list)):
    
    l_w = lw_list[l_w_idx]
    l_c = 1
    
    # go to daughter folder
    folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
    try:
        os.chdir(mother_folder + "\\"+ folder_name) #windows
    except:
        os.chdir(mother_folder + "/"+ folder_name) #linux
    
    time_points = np.loadtxt("time_points.txt", delimiter=',')
    
    fw_beta_a_w_integ = np.loadtxt("fw_beta_a_w_integ.txt", delimiter=',')
    Wa_beta_a_w_integ = np.loadtxt("Wa_beta_a_w_integ.txt", delimiter=',')
    w_a_mat = np.loadtxt("w_a_mat.txt", delimiter=',')
    w_mat = np.loadtxt("w_mat.txt", delimiter=',')
    fw_mat = w_a_mat/w_mat
    org_names = []
    with open("org_names.txt", "r") as f:
        for line in f:
            org_names.append(line.strip())
    
    # eliminating outliers
    outlier_list = []
    for outlier_name in outlier_names_list:
        if outlier_name in org_names:
            outlier_list.append(  org_names.index(outlier_name)  )
    fw_beta_a_w_integ = np.delete(fw_beta_a_w_integ, outlier_list, axis=0)
    Wa_beta_a_w_integ = np.delete(Wa_beta_a_w_integ, outlier_list, axis=0)
    w_a_mat = np.delete(w_a_mat, outlier_list, axis=0)
    w_mat = np.delete(w_mat, outlier_list, axis=0)
    fw_mat = np.delete(fw_mat, outlier_list, axis=0)
    # eliminating outliers
    
    n_orgs = np.shape(fw_beta_a_w_integ)[0]
    
    f_beta_w_integ_avg_mat[l_w_idx,:] = np.mean(fw_beta_a_w_integ, axis=0)
    f_beta_w_integ_err_mat[l_w_idx,:] = np.std(fw_beta_a_w_integ, axis=0)/np.sqrt(n_orgs-1)
    N_beta_w_integ_avg_mat[l_w_idx,:] = np.mean(Wa_beta_a_w_integ, axis=0)
    N_beta_w_integ_err_mat[l_w_idx,:] = np.std(Wa_beta_a_w_integ, axis=0)/np.sqrt(n_orgs-1)
    
    w_a_avg_all_lw[l_w_idx,:] = np.mean(w_a_mat, axis=0)
    w_a_err_all_lw[l_w_idx,:] = np.std(w_a_mat, axis=0)/np.sqrt(n_orgs-1)
    fw_avg_all_lw[l_w_idx,:] = np.mean(fw_mat, axis=0)
    fw_err_all_lw[l_w_idx,:] = np.std(fw_mat, axis=0)/np.sqrt(n_orgs-1)
    
    # fitting by f
    x = time_points
    avg = f_beta_w_integ_avg_mat[l_w_idx,:]
    err = f_beta_w_integ_err_mat[l_w_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    F_w_avg[l_w_idx,:] = poly(fit_plot_x)
    dF_dx = np.gradient(F_w_avg[l_w_idx,:], fit_plot_x)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, F_w_avg[l_w_idx,:])
    plt.show()
    
    avg = fw_avg_all_lw[l_w_idx,:]
    err = fw_err_all_lw[l_w_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    weights[err==0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=1, w=weights)
    poly = np.poly1d(coeffs)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.show()
    
    beta_w_a_avg_by_f[l_w_idx,:] = dF_dx / poly(fit_plot_x)
    # beta_w_a_err_by_f[l_w_idx,:] = 
    # fitting by f
    
    # fitting by N
    x = time_points
    avg = N_beta_w_integ_avg_mat[l_w_idx,:]
    err = N_beta_w_integ_err_mat[l_w_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    G_w_avg[l_w_idx,:] = poly(fit_plot_x)
    dG_dx = np.gradient(G_w_avg[l_w_idx,:], fit_plot_x)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.show()
    
    avg = w_a_avg_all_lw[l_w_idx,:]
    err = w_a_err_all_lw[l_w_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=1, w=weights)
    poly = np.poly1d(coeffs)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.show()
    
    beta_w_a_avg_by_N[l_w_idx,:] = dG_dx / poly(fit_plot_x)
    # beta_w_a_err_by_f[l_w_idx,:] = 
    # fitting by N
    
    os.chdir(mother_folder)
# W data reading #

# W plot and fit #
plt.figure()
x = time_points
for l_w_idx in range(len(lw_list)):
    l_w = lw_list[l_w_idx]
    
    avg = f_beta_w_integ_avg_mat[l_w_idx,:]
    err = f_beta_w_integ_err_mat[l_w_idx,:]
    # plt.plot(x, avg, label='lw='+str(int(l_w)))
    plt.errorbar(x, avg, yerr=err , label='lw='+str(int(l_w)), linestyle=None, fmt='o')
        
    # weights = 1/err
    # weights[0] = 1/(1e-4)
    # coeffs = np.polyfit(x, avg, deg=3, w=weights)
    # poly = np.poly1d(coeffs)
    # fit_plot_x = np.linspace(np.min(x), np.max(x), 1000)
    # plt.plot(fit_plot_x, poly(fit_plot_x), linestyle='--', color='k')
    plt.plot(fit_plot_x, F_w_avg[l_w_idx,:], linestyle='--', color='k')
    # dy_dx = np.gradient(poly(fit_plot_x), fit_plot_x)
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\int_0^t\langle\beta^a_W(s) f(s)\rangle ds$', fontsize=20)
plt.tight_layout()
plt.savefig("f_beta_w"+".png", dpi=300)
# plt.close()

plt.figure()
x = time_points
for l_w_idx in range(len(lw_list)):
    l_w = lw_list[l_w_idx]
    avg = N_beta_w_integ_avg_mat[l_w_idx,:]
    err = N_beta_w_integ_err_mat[l_w_idx,:]
    # plt.plot(x, avg, label='lw='+str(int(l_w)))
    plt.errorbar(x, avg, yerr=err , label='lw='+str(int(l_w)), linestyle=None, fmt='o')
    # weights = 1/err
    # weights[0] = 1/(1e-4)
    # coeffs = np.polyfit(x, avg, deg=3, w=weights)
    # poly = np.poly1d(coeffs)
    # fit_plot_x = np.linspace(np.min(x), np.max(x), 1000)
    # plt.plot(fit_plot_x, poly(fit_plot_x), linestyle='--', color='k')
    plt.plot(fit_plot_x, G_w_avg[l_w_idx,:], linestyle='--', color='k')
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\int_0^t\langle\beta^a_W(s) W^a(s)\rangle ds$', fontsize=20)
plt.grid()
plt.tight_layout()
plt.savefig("N_beta_w"+".png", dpi=300)
# plt.close()

plt.figure()
x = time_points
colors = [
    'blue', 'orange', 'green', 'red', 'purple',
    'brown', 'pink', 'gray', 'olive', 'cyan'
]
for l_w_idx in range(len(lw_list)):
    l_w = lw_list[l_w_idx]
    
    plt.plot(fit_plot_x, beta_w_a_avg_by_f[l_w_idx,:], label='lw='+str(int(l_w))+', f', linestyle='--', color=colors[l_w_idx])
    plt.plot(fit_plot_x, beta_w_a_avg_by_N[l_w_idx,:], label='lw='+str(int(l_w))+', N', color=colors[l_w_idx])

plt.plot(fit_plot_x, 0*fit_plot_x+beta_w_0, label=r'$\beta^0_W$')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\langle\beta^a_W(t)\rangle$', fontsize=20)
plt.grid()
plt.tight_layout()
plt.savefig("beta_w"+".png", dpi=300)
# plt.close()
# W plot and fit #




# C data reading #
for l_c_idx in range(len(lc_list)):
    
    l_c = lc_list[l_c_idx]
    l_w = 1
    
    # go to daughter folder
    folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
    try:
        os.chdir(mother_folder + "\\"+ folder_name) #windows
    except:
        os.chdir(mother_folder + "/"+ folder_name) #linux
    
    time_points = np.loadtxt("time_points.txt", delimiter=',')
    
    
    fc_beta_a_c_integ = np.loadtxt("fc_beta_a_c_integ.txt", delimiter=',')
    Ca_beta_a_c_integ = np.loadtxt("Ca_beta_a_c_integ.txt", delimiter=',')
    c_a_mat = np.loadtxt("c_a_mat.txt", delimiter=',')
    c_mat = np.loadtxt("c_mat.txt", delimiter=',')
    fc_mat = c_a_mat/c_mat
    org_names = []
    with open("org_names.txt", "r") as f:
        for line in f:
            org_names.append(line.strip())
    
    # eliminating outliers
    outlier_list = []
    for outlier_name in outlier_names_list:
        if outlier_name in org_names:
            outlier_list.append(  org_names.index(outlier_name)  )
    fc_beta_a_c_integ = np.delete(fc_beta_a_c_integ, outlier_list, axis=0)
    Ca_beta_a_c_integ = np.delete(Ca_beta_a_c_integ, outlier_list, axis=0)
    c_a_mat = np.delete(c_a_mat, outlier_list, axis=0)
    c_mat = np.delete(c_mat, outlier_list, axis=0)
    fc_mat = np.delete(fc_mat, outlier_list, axis=0)
    # eliminating outliers
    
    n_orgs = np.shape(fc_beta_a_c_integ)[0]
    
    f_beta_c_integ_avg_mat[l_c_idx,:] = np.mean(fc_beta_a_c_integ, axis=0)
    f_beta_c_integ_err_mat[l_c_idx,:] = np.std(fc_beta_a_c_integ, axis=0)/np.sqrt(n_orgs-1)
    N_beta_c_integ_avg_mat[l_c_idx,:] = np.mean(Ca_beta_a_c_integ, axis=0)
    N_beta_c_integ_err_mat[l_c_idx,:] = np.std(Ca_beta_a_c_integ, axis=0)/np.sqrt(n_orgs-1)
    
    c_a_avg_all_lc[l_c_idx,:] = np.mean(c_a_mat, axis=0)
    c_a_err_all_lc[l_c_idx,:] = np.std(c_a_mat, axis=0)/np.sqrt(n_orgs-1)
    fc_avg_all_lc[l_c_idx,:] = np.mean(fc_mat, axis=0)
    fc_err_all_lc[l_c_idx,:] = np.std(fc_mat, axis=0)/np.sqrt(n_orgs-1)
    
    # fitting by f
    x = time_points
    avg = f_beta_c_integ_avg_mat[l_c_idx,:]
    err = f_beta_c_integ_err_mat[l_c_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    F_c_avg[l_c_idx,:] = poly(fit_plot_x)
    dF_dx = np.gradient(F_c_avg[l_c_idx,:], fit_plot_x)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, F_c_avg[l_c_idx,:])
    plt.title('int_fc_b_a_c')
    plt.show()
    
    avg = fc_avg_all_lc[l_c_idx,:]
    err = fc_err_all_lc[l_c_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    weights[err==0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.title('fc')
    plt.show()
    
    beta_c_a_avg_by_f[l_c_idx,:] = dF_dx / poly(fit_plot_x)
    # fitting by f
    
    # fitting by N
    x = time_points
    avg = N_beta_c_integ_avg_mat[l_c_idx,:]
    err = N_beta_c_integ_err_mat[l_c_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    G_c_avg[l_c_idx,:] = poly(fit_plot_x)
    dG_dx = np.gradient(G_c_avg[l_c_idx,:], fit_plot_x)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.title('int_N_b_a_c')
    plt.show()
    
    
    avg = c_a_avg_all_lc[l_c_idx,:]
    err = c_a_err_all_lc[l_c_idx,:]
    weights = 1/err
    weights[0] = 1/(1e-4)
    coeffs = np.polyfit(x, avg, deg=3, w=weights)
    poly = np.poly1d(coeffs)
    
    plt.figure()
    plt.errorbar(x, avg, yerr=err ,  linestyle=None, fmt='o')
    plt.plot(fit_plot_x, poly(fit_plot_x))
    plt.yscale('log')
    plt.title('c_a')
    plt.show()
    
    # plt.figure()
    # # plt.plot(fit_plot_x, dG_dx-beta_c_0*poly(fit_plot_x))
    # # plt.plot(fit_plot_x, dG_dx)
    # plt.plot(fit_plot_x, poly(fit_plot_x))
    # plt.title('svdsvsf')
    # plt.show()
    
    beta_c_a_avg_by_N[l_c_idx,:] = dG_dx / poly(fit_plot_x)
    # beta_w_a_err_by_f[l_w_idx,:] = 
    # fitting by N
    
    os.chdir(mother_folder)
# C data reading #

# C plot and fit #
plt.figure()
x = time_points
for l_c_idx in range(len(lc_list)):
    l_c = lc_list[l_c_idx]
    
    avg = f_beta_c_integ_avg_mat[l_c_idx,:]
    err = f_beta_c_integ_err_mat[l_c_idx,:]
    # plt.plot(x, avg, label='lw='+str(int(l_w)))
    plt.errorbar(x, avg, yerr=err , label='lc='+str(int(l_c)), linestyle=None, fmt='o')
        
    # weights = 1/err
    # weights[0] = 1/(1e-4)
    # coeffs = np.polyfit(x, avg, deg=3, w=weights)
    # poly = np.poly1d(coeffs)
    # fit_plot_x = np.linspace(np.min(x), np.max(x), 1000)
    
    plt.plot(fit_plot_x, F_c_avg[l_c_idx,:], linestyle='--', color='k')
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\int_0^t\langle\beta^a_C(s) f(s)\rangle ds$', fontsize=20)

plt.tight_layout()
plt.savefig("f_beta_C"+".png", dpi=300)
# plt.close()

plt.figure()
x = time_points
for l_c_idx in range(len(lc_list)):
    l_c = lc_list[l_c_idx]
    
    avg = N_beta_c_integ_avg_mat[l_c_idx,:]
    err = N_beta_c_integ_err_mat[l_c_idx,:]
    # plt.plot(x, avg, label='lw='+str(int(l_w)))
    plt.errorbar(x, avg, yerr=err , label='lc='+str(int(l_c)), linestyle=None, fmt='o')
        
    # weights = 1/err
    # weights[0] = 1/(1e-4)
    # coeffs = np.polyfit(x, avg, deg=3, w=weights)
    # poly = np.poly1d(coeffs)
    # fit_plot_x = np.linspace(np.min(x), np.max(x), 1000)
    
    plt.plot(fit_plot_x, G_c_avg[l_c_idx,:], linestyle='--', color='k')
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\int_0^t\langle\beta^a_C(s) C^a(s)\rangle ds$', fontsize=20)
plt.grid()
plt.tight_layout()
plt.savefig("N_beta_c"+".png", dpi=300)
# plt.close()

plt.figure()
x = time_points
colors = [
    'blue', 'orange', 'green', 'red', 'purple',
    'brown', 'pink', 'gray', 'olive', 'cyan'
]
for l_c_idx in range(len(lc_list)):
    l_c = lc_list[l_c_idx]
    
    plt.plot(fit_plot_x, beta_c_a_avg_by_f[l_c_idx,:], label='lc='+str(int(l_c))+', f', linestyle='--', color=colors[l_c_idx])
    plt.plot(fit_plot_x, beta_c_a_avg_by_N[l_c_idx,:], label='lc='+str(int(l_c))+', N', color=colors[l_c_idx])

plt.plot(fit_plot_x, 0*fit_plot_x+beta_c_0, label=r'$\beta^0_C$')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('time (h)', fontsize=15)
plt.ylabel(r'$\langle\beta^a_C(t)\rangle$', fontsize=20)
plt.grid()
plt.tight_layout()
plt.savefig("beta_c"+".png", dpi=300)
# plt.close()
# C plot and fit #