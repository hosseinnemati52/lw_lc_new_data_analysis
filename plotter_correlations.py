# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:21:47 2026

@author: Nemat002
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import os
from scipy.stats import pearsonr


def time_plotter_func(x, y_mat, ylabel, file_name, color, switch_all_plot):
    
    if np.std(y_mat[:,0])<0.001 and (  abs(np.mean(y_mat[:,0])-1)<0.001 ):
        switch_norm=1
    else:
        switch_norm=0
        
    plt.figure()
    if switch_all_plot:
        for i in range(n_orgs):
            plt.plot(x, y_mat[i,:], color=color)
    plt.xlabel('time (h)', fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if switch_norm:
        plt.yscale("log")
    avg = np.mean(y_mat, axis=0)
    err = np.std(y_mat, axis=0)/np.sqrt(n_orgs-1)
    plt.plot(x, avg, color='k')
    plt.errorbar(x, avg, yerr=err, color='k')
    
    
    if file_name=='fw_beta_a_w_integ' or file_name=='Wa_beta_a_w_integ':
        
        weights = 1/err
        weights[0] = 1/(1e-4)
        coeffs = np.polyfit(x, avg, deg=3, w=weights)
        
        poly = np.poly1d(coeffs)

        fit_plot_x = np.linspace(np.min(x), np.max(x), 1000)
        plt.plot(fit_plot_x, poly(fit_plot_x), linestyle='--', color='brown')
        

    
    plt.tight_layout()
    plt.savefig(file_name+".png", dpi=300)
    plt.close()
    
    return 0

def scatter_plotter(x, y, xlabel, ylabel, file_name, color):
    
    corr_mat = np.corrcoef(x,y)
    
    r, p = pearsonr(x, y)
    
    plt.figure()
    plt.scatter(x, y, color=color)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    title = 'corr: '+str( round(r, 3) )+', p_val: '+str(round(p, 5))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(file_name+".png", dpi=300)
    plt.close()
    
    
    return r, p

def WLS_fitter(x_fit, y_fit, yerr_fit):
    
    ######################################
    # fitting with uncertainty
    x_lin_fit =  x_fit
    y_lin_fit =  y_fit
    y_lin_ERR_fit =  yerr_fit
    
    Delta_mat =     np.array([[0., 0.], [0., 0.]])
    intercept_mat = np.array([[0., 0.], [0., 0.]])
    slope_mat =     np.array([[0., 0.], [0., 0.]])
    
    Delta_mat[0,0] = np.sum(1/y_lin_ERR_fit**2)
    Delta_mat[0,1] = np.sum(x_lin_fit/y_lin_ERR_fit**2)
    Delta_mat[1,0] = np.sum(x_lin_fit/y_lin_ERR_fit**2)
    Delta_mat[1,1] = np.sum(x_lin_fit**2/y_lin_ERR_fit**2)
    Delta = np.linalg.det(Delta_mat)
    
    intercept_mat[0,0] = np.sum(y_lin_fit/y_lin_ERR_fit**2)
    intercept_mat[0,1] = np.sum(x_lin_fit/y_lin_ERR_fit**2)
    intercept_mat[1,0] = np.sum(x_lin_fit*y_lin_fit/y_lin_ERR_fit**2)
    intercept_mat[1,1] = np.sum(x_lin_fit**2/y_lin_ERR_fit**2)
    intercept = (1/Delta)*np.linalg.det(intercept_mat)
    
    slope_mat[0,0] = np.sum(1/y_lin_ERR_fit**2)
    slope_mat[0,1] = np.sum(y_lin_fit/y_lin_ERR_fit**2)
    slope_mat[1,0] = np.sum(x_lin_fit/y_lin_ERR_fit**2)
    slope_mat[1,1] = np.sum(x_lin_fit*y_lin_fit/y_lin_ERR_fit**2)
    slope = (1/Delta)*np.linalg.det(slope_mat)
    
    
    intercept_ERR = np.sqrt(  (1/Delta) * np.sum(x_lin_fit**2/y_lin_ERR_fit**2) )
    slope_ERR =     np.sqrt(  (1/Delta) * np.sum(1.0         /y_lin_ERR_fit**2) )
    ######################################
    
    return slope, slope_ERR, intercept, intercept_ERR

lw_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]

# lw_list = [20]
# lc_list = [20]

# data for 65 and 70
WT_pure = np.genfromtxt("WT_pure.csv", delimiter=",", filling_values=np.nan)[:,1:]
C_pure = np.genfromtxt("C_pure.csv", delimiter=",", filling_values=np.nan)[:,1:]
WT_mix = np.genfromtxt("WT_mix.csv", delimiter=",", filling_values=np.nan)[:,1:]
C_mix = np.genfromtxt("C_mix.csv", delimiter=",", filling_values=np.nan)[:,1:]

log_WT_mix_norm = np.log(WT_mix / WT_mix[[0],:])
log_C_mix_norm = np.log(C_mix / C_mix[[0],:])
delta_WT_mix = WT_mix - WT_mix[[0],:]
delta_C_mix = C_mix - C_mix[[0],:]

#
log_w_to_w0 = dict()
all_data = log_WT_mix_norm[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
log_w_to_w0[70] = np.array([mean_val, SEM_val])

all_data = log_WT_mix_norm[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
log_w_to_w0[65] = np.array([mean_val, SEM_val])

log_c_to_c0 = dict()
all_data = log_C_mix_norm[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
log_c_to_c0[70] = np.array([mean_val, SEM_val])

all_data = log_C_mix_norm[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
log_c_to_c0[65] = np.array([mean_val, SEM_val])
#

#
w_minus_w0 = dict()
all_data = delta_WT_mix[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
w_minus_w0[70] = np.array([mean_val, SEM_val])

all_data = delta_WT_mix[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
w_minus_w0[65] = np.array([mean_val, SEM_val])


c_minus_c0 = dict()
all_data = delta_C_mix[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
c_minus_c0[70] = np.array([mean_val, SEM_val])

all_data = delta_C_mix[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
c_minus_c0[65] = np.array([mean_val, SEM_val])
#

#
w_extra = dict()
all_data = WT_mix[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
w_extra[70] = np.array([mean_val, SEM_val])

all_data = WT_mix[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
w_extra[65] = np.array([mean_val, SEM_val])

c_extra = dict()
all_data = C_mix[-1,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
c_extra[70] = np.array([mean_val, SEM_val])

all_data = C_mix[-2,:]
all_data = all_data[~np.isnan(all_data)]
mean_val = np.mean(all_data)
SEM_val = np.std(all_data)/np.sqrt(len(all_data)-1)
c_extra[65] = np.array([mean_val, SEM_val])
#



mother_folder = os.getcwd()
os.chdir(mother_folder)

pure_betas = np.loadtxt('pure_organoid_exponents.csv', delimiter=',')
beta_w_0 = pure_betas[0,0]
beta_c_0 = pure_betas[0,1]

conc_thresh = 150

for l_w in lw_list:
    for l_c in lc_list:
        
        print('lw: '+str(l_w)+', '+'lc: '+str(l_c))
        
        # go to daughter folder
        folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        try:
            os.chdir(mother_folder + "\\"+ folder_name) #windows
        except:
            os.chdir(mother_folder + "/"+ folder_name) #linux
        
        org_names = []
        del org_names
        org_names = []
        with open("org_names.txt", "r") as f:
            for line in f:
                org_names.append(line.strip())
        
        # reading data
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
        
        # eliminating outliers
        outlier_names_list = ['E-139_S-18_O-1']
        outlier_list = []
        for outlier_name in outlier_names_list:
            if outlier_name in org_names:
                outlier_list.append(  org_names.index(outlier_name)  )
        w_mat = np.delete(w_mat, outlier_list, axis=0)
        w_b_mat = np.delete(w_b_mat, outlier_list, axis=0)
        w_a_mat = np.delete(w_a_mat, outlier_list, axis=0)
        w_v_mat = np.delete(w_v_mat, outlier_list, axis=0)
        w_u_mat = np.delete(w_u_mat, outlier_list, axis=0)
        
        c_mat = np.delete(c_mat, outlier_list, axis=0)
        c_b_mat = np.delete(c_b_mat, outlier_list, axis=0)
        c_a_mat = np.delete(c_a_mat, outlier_list, axis=0)
        c_v_mat = np.delete(c_v_mat, outlier_list, axis=0)
        c_u_mat = np.delete(c_u_mat, outlier_list, axis=0)
        
        avg_num_c_visib_to_aff_w_mat = np.delete(avg_num_c_visib_to_aff_w_mat, outlier_list, axis=0)
        avg_num_w_visib_to_aff_c_mat = np.delete(avg_num_w_visib_to_aff_c_mat, outlier_list, axis=0)
        
        vol_lum = np.delete(vol_lum, outlier_list, axis=0)
        
        fw_beta_a_w_integ = np.delete(fw_beta_a_w_integ, outlier_list, axis=0)
        fc_beta_a_c_integ = np.delete(fc_beta_a_c_integ, outlier_list, axis=0)
        
        Wa_beta_a_w_integ = np.delete(Wa_beta_a_w_integ, outlier_list, axis=0)
        Ca_beta_a_c_integ = np.delete(Ca_beta_a_c_integ, outlier_list, axis=0)
        
        # eliminating outliers
        
        fw_mat = w_a_mat/w_mat # fraction of affected
        fc_mat = c_a_mat/c_mat
        
        dt = time_points[1] - time_points[0]
        
        n_orgs = np.shape(w_mat)[0]
        n_time = np.shape(w_mat)[1]
        
        # beta_w_a = np.loadtxt("beta_w_a.txt", delimiter=',')
        # beta_c_a = np.loadtxt("beta_c_a.txt", delimiter=',')
        
        # temporal plots
        w_color = 'm'
        c_color = 'g'
        
        fff =1
        if fff:
            
            time_plotter_func(time_points, w_mat,   'w',   'w_tot', w_color,1)
            time_plotter_func(time_points, w_b_mat, 'w_b', 'w_b',   w_color,1)
            time_plotter_func(time_points, w_a_mat, 'w_a', 'w_a',   w_color,1)
            time_plotter_func(time_points, w_v_mat, 'w_v', 'w_v',   w_color,1)
            time_plotter_func(time_points, w_u_mat, 'w_u', 'w_u',   w_color,1)
            time_plotter_func(time_points, avg_num_c_visib_to_aff_w_mat, 'avg_num_c_visib_to_aff_w_mat', 'avg_num_c_visib_to_aff_w_mat',   w_color,1)
            time_plotter_func(time_points, fw_mat, 'fw_mat', 'fw_mat',   w_color,1)
            
            time_plotter_func(time_points, w_mat/w_mat[:, [0]],   'w_norm',   'w_tot_norm', w_color,1)
            time_plotter_func(time_points, w_b_mat/w_b_mat[:, [0]], 'w_b_norm', 'w_b_norm',   w_color,1)
            time_plotter_func(time_points, w_a_mat/w_a_mat[:, [0]], 'w_a_norm', 'w_a_norm',   w_color,1)
            time_plotter_func(time_points, w_v_mat/w_v_mat[:, [0]], 'w_v_norm', 'w_v_norm',   w_color,1)
            time_plotter_func(time_points, w_u_mat/w_u_mat[:, [0]], 'w_u_norm', 'w_u_norm',   w_color,1)
            time_plotter_func(time_points, avg_num_c_visib_to_aff_w_mat/avg_num_c_visib_to_aff_w_mat[:, [0]], 'avg_num_c_visib_to_aff_w_mat_norm', 'avg_num_c_visib_to_aff_w_mat_norm',   w_color,1)
            
            time_plotter_func(time_points, c_mat,   'c',   'c_tot', c_color,1)
            time_plotter_func(time_points, c_b_mat, 'c_b', 'c_b',   c_color,1)
            time_plotter_func(time_points, c_a_mat, 'c_a', 'c_a',   c_color,1)
            time_plotter_func(time_points, c_v_mat, 'c_v', 'c_v',   c_color,1)
            time_plotter_func(time_points, c_u_mat, 'c_u', 'c_u',   c_color,1)
            time_plotter_func(time_points, avg_num_w_visib_to_aff_c_mat, 'avg_num_w_visib_to_aff_c_mat', 'avg_num_w_visib_to_aff_c_mat',   c_color,1)
            time_plotter_func(time_points, fc_mat, 'fc_mat', 'fc_mat',   c_color,1)
            
            time_plotter_func(time_points, c_mat/c_mat[:, [0]],   'c_norm',   'c_tot_norm', c_color,1)
            time_plotter_func(time_points, c_b_mat/c_b_mat[:, [0]], 'c_b_norm', 'c_b_norm',   c_color,1)
            time_plotter_func(time_points, c_a_mat/c_a_mat[:, [0]], 'c_a_norm', 'c_a_norm',   c_color,1)
            time_plotter_func(time_points, c_v_mat/c_v_mat[:, [0]], 'c_v_norm', 'c_v_norm',   c_color,1)
            time_plotter_func(time_points, c_u_mat/c_u_mat[:, [0]], 'c_u_norm', 'c_u_norm',   c_color,1)
            time_plotter_func(time_points, avg_num_w_visib_to_aff_c_mat/avg_num_w_visib_to_aff_c_mat[:, [0]], 'avg_num_w_visib_to_aff_c_mat_norm', 'avg_num_w_visib_to_aff_c_mat_norm',   c_color,1)
            
            time_plotter_func(time_points, vol_lum, 'vol_lum', 'vol_lum',   'b',1)
            time_plotter_func(time_points, vol_lum/vol_lum[:, [0]], 'vol_lum_norm', 'vol_lum_norm',   'b',1)
            
            # time_plotter_func(time_points, fw_beta_a_w, 'fw_beta_a_w', 'fw_beta_a_w',   w_color)
            # time_plotter_func(time_points, fc_beta_a_c, 'fc_beta_a_c', 'fc_beta_a_c',   c_color)
            
            time_plotter_func(time_points, fw_beta_a_w_integ, 'fw_beta_a_w_integ', 'fw_beta_a_w_integ',   w_color,0)
            time_plotter_func(time_points, fc_beta_a_c_integ, 'fc_beta_a_c_integ', 'fc_beta_a_c_integ',   c_color,0)
            
            time_plotter_func(time_points, Wa_beta_a_w_integ, 'Wa_beta_a_w_integ', 'Wa_beta_a_w_integ',   w_color,0)
            time_plotter_func(time_points, Ca_beta_a_c_integ, 'Ca_beta_a_c_integ', 'Ca_beta_a_c_integ',   c_color,0)
            
            # # plotting with extra data points
            # len_extra = n_time+2
            # fw_beta_a_w_integ_extra = np.zeros(len_extra)
            # fc_beta_a_c_integ_extra = np.zeros(len_extra)
            # Wa_beta_a_w_integ_extra = np.zeros(len_extra)
            # Ca_beta_a_c_integ_extra = np.zeros(len_extra)
            
            # fw_beta_a_w_integ_extra_err = np.zeros(len_extra)
            # fc_beta_a_c_integ_extra_err = np.zeros(len_extra)
            # Wa_beta_a_w_integ_extra_err = np.zeros(len_extra)
            # Ca_beta_a_c_integ_extra_err = np.zeros(len_extra)
            
            # fw_beta_a_w_integ_extra[:n_time] = np.mean(fw_beta_a_w_integ, axis=0)
            # fc_beta_a_c_integ_extra[:n_time] = np.mean(fc_beta_a_c_integ, axis=0)
            # Wa_beta_a_w_integ_extra[:n_time] = np.mean(Wa_beta_a_w_integ, axis=0)
            # Ca_beta_a_c_integ_extra[:n_time] = np.mean(Ca_beta_a_c_integ, axis=0)
            
            # fw_beta_a_w_integ_extra_err[:n_time] = np.std(fw_beta_a_w_integ, axis=0)/np.sqrt(n_orgs-1)
            # fc_beta_a_c_integ_extra_err[:n_time] = np.std(fc_beta_a_c_integ, axis=0)/np.sqrt(n_orgs-1)
            # Wa_beta_a_w_integ_extra_err[:n_time] = np.std(Wa_beta_a_w_integ, axis=0)/np.sqrt(n_orgs-1)
            # Ca_beta_a_c_integ_extra_err[:n_time] = np.std(Ca_beta_a_c_integ, axis=0)/np.sqrt(n_orgs-1)
            
            
            # # W, C
            # n_fit_data = 3
            # data_for_fit = np.mean(w_a_mat, axis=0)[-n_fit_data:]
            # err_for_fit  = np.std(w_a_mat, axis=0)[-n_fit_data:]/np.sqrt(n_orgs-1)
            # slope, slope_ERR, intercept, intercept_ERR = WLS_fitter(time_points[-n_fit_data:], data_for_fit, err_for_fit)
            # W_a_avg_fitted = np.array([65, 70])* slope + intercept
            
            # data_for_fit = np.std(w_a_mat, axis=0)[-n_fit_data:]/np.sqrt(n_orgs-1)
            # err_for_fit  = 0.0 * data_for_fit + (1e-7)
            # slope, slope_ERR, intercept, intercept_ERR = WLS_fitter(time_points[-n_fit_data:], data_for_fit, err_for_fit)
            # W_a_err_fitted = np.array([65, 70])* slope + intercept
            
            # W_u_avg_fitted = np.array([w_extra[65][0], w_extra[70][0]]) - W_a_avg_fitted
            # W_u_err_fitted = np.array([w_extra[65][1], w_extra[70][1]]) + W_a_err_fitted
            
            
            # # W, C
            
            # # f
            # # f
            
            # ggg
            # # plotting with extra data points
            
            # temporal plots
        
        # rar
        # scatter plots
        corr_list = []
        p_val_list = []
        t_eval = 60
        # y = beta_w_0*t_eval - np.log(w_mat[:,-1] / w_mat[:,0])
        y_w_a = w_mat[:,-1] - w_mat[:,0] - beta_w_0 * dt * (   np.sum(w_u_mat[:,1:-1], axis=1)  + 0.5*(w_u_mat[:,0]+w_u_mat[:,-1])   )
        x_w_a =                            beta_w_0 * dt * (   np.sum(w_a_mat[:,1:-1], axis=1)  + 0.5*(w_a_mat[:,0]+w_a_mat[:,-1])   )
        y = y_w_a / x_w_a
        
        x = dt * ( np.sum(c_v_mat[:,1:-1], axis=1) + 0.5*(c_v_mat[:,0]+c_v_mat[:,-1]) )
        # corr = scatter_plotter(x, y, r'$\int_0^{60}C^v(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_Cv', w_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60}C^v(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_Cv', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        x = np.sum(c_v_mat[:,1:-1]/w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(c_v_mat[:,0]/w_a_mat[:,0]+c_v_mat[:,-1]/w_a_mat[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60}C^v(s)/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_CvWa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60}C^v(s)/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_CvWa', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        # x = np.sum(1/w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(1/w_a_mat[:,0]+1/w_a_mat[:,-1])
        x = np.sum(w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(w_a_mat[:,0]+w_a_mat[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_Wa', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        # x = np.sum(1/w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(1/w_a_mat[:,0]+1/w_a_mat[:,-1])
        x = c_v_mat[:,0]
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$C^v(0)$' , r'$y_a^w/x_a^w$', 'cor_W_Cv0', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        # x = np.sum(1/w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(1/w_a_mat[:,0]+1/w_a_mat[:,-1])
        x = c_v_mat[:,0]/w_a_mat[:,0]
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$C^v(0)/W^a(0)$' , r'$y_a^w/x_a^w$', 'cor_W_Cv0Wa0', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        integrand = avg_num_c_visib_to_aff_w_mat
        x = np.sum(integrand[:,1:-1], axis=1)*dt + 0.5*dt*(integrand[:,0]+integrand[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} n_d(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_nd', w_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        if l_w > conc_thresh:
            y = w_mat[:,-1] / w_mat[:,0]
            cumul_c_v_mat = 0.0 * c_v_mat.copy()
            for org_c in range(n_orgs):
                for t_c in range(1, n_time):
                    cumul_c_v_mat[org_c, t_c] = 0.5*dt*c_v_mat[org_c, t_c] + dt*np.sum(c_v_mat[org_c, 1:t_c]) + 0.5*dt*c_v_mat[org_c, 0]
            integrand = cumul_c_v_mat/vol_lum
            x = np.sum(integrand[:,1:-1], axis=1)*dt + 0.5*dt*(integrand[:,0]+integrand[:,-1])
            # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
            corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} CC^v(s)/V_{\ell}(s) ds$' , r'$W(60)/W(0)$', 'cor_W_conc', w_color)
            corr_list.append(corr)
            p_val_list.append(p_val)
        else:
            corr_list.append(0.0)
            p_val_list.append(0.0)
            
        
        
        # y = np.log(c_mat[:,-1] / c_mat[:,0]) - beta_c_0*t_eval
        y_c_a = c_mat[:,-1] - c_mat[:,0] - beta_c_0 * dt * (   np.sum(c_u_mat[:,1:-1], axis=1)  + 0.5*(c_u_mat[:,0]+c_u_mat[:,-1])   )
        x_c_a =                            beta_c_0 * dt * (   np.sum(c_a_mat[:,1:-1], axis=1)  + 0.5*(c_a_mat[:,0]+c_a_mat[:,-1])   )
        y = y_c_a / x_c_a
        
        x = np.sum(w_v_mat[:,1:-1], axis=1)*dt + 0.5*dt*(w_v_mat[:,0]+w_v_mat[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60}W^v(s) ds$' , r'$-\beta_0^c t + \log\frac{C(60)}{C(0)}$', 'cor_C_Wv', c_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60}W^v(s) ds$' , r'$y_a^c/x_a^c$', 'cor_C_Wv', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        x = np.sum(w_v_mat[:,1:-1]/c_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(w_v_mat[:,0]/c_a_mat[:,0]+w_v_mat[:,-1]/c_a_mat[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60}W^v(s)/C^a(s) ds$' , r'$-\beta_0^c t + \log\frac{C(60)}{C(0)}$', 'cor_C_WvCa', c_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60}W^v(s)/C^a(s) ds$' , r'$y_a^c/x_a^c$', 'cor_C_WvCa', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        # x = np.sum(1/c_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(1/c_a_mat[:,0]+1/c_a_mat[:,-1])
        x = np.sum(c_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(c_a_mat[:,0]+c_a_mat[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/C^a(s) ds$' , r'$-\beta_0^c t + \log\frac{C(60)}{C(0)}$', 'cor_C_1Ca', c_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/C^a(s) ds$' , r'$y_a^c/x_a^c$', 'cor_C_1Ca', c_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} C^a(s) ds$' , r'$y_a^c/x_a^c$', 'cor_C_Ca', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        x = w_v_mat[:,0]
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$W^v(0)$' , r'$y_a^c/x_a^c$', 'cor_C_Wv0', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        # x = np.sum(1/w_a_mat[:,1:-1], axis=1)*dt + 0.5*dt*(1/w_a_mat[:,0]+1/w_a_mat[:,-1])
        x = w_v_mat[:,0]/c_a_mat[:,0]
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$y_a^w/x_a^w$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$W^v(0)/C^a(0)$' , r'$y_a^c/x_a^c$', 'cor_C_Wv0Ca0', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        integrand = avg_num_w_visib_to_aff_c_mat
        x = np.sum(integrand[:,1:-1], axis=1)*dt + 0.5*dt*(integrand[:,0]+integrand[:,-1])
        # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
        corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} n_d(s) ds$' , r'$y_a^c/x_a^c$', 'cor_C_nd', c_color)
        corr_list.append(corr)
        p_val_list.append(p_val)
        
        if l_c > conc_thresh:
            y = c_mat[:,-1] / c_mat[:,0]
            cumul_w_v_mat = 0.0 * w_v_mat.copy()
            for org_c in range(n_orgs):
                for t_c in range(1, n_time):
                    cumul_w_v_mat[org_c, t_c] = 0.5*dt*w_v_mat[org_c, t_c] + dt*np.sum(w_v_mat[org_c, 1:t_c]) + 0.5*dt*w_v_mat[org_c, 0]
            integrand = cumul_w_v_mat/vol_lum
            x = np.sum(integrand[:,1:-1], axis=1)*dt + 0.5*dt*(integrand[:,0]+integrand[:,-1])
            # corr = scatter_plotter(x, y, r'$\int_0^{60} 1/W^a(s) ds$' , r'$\beta_0^w t - \log\frac{W(60)}{W(0)}$', 'cor_W_1Wa', w_color)
            corr, p_val = scatter_plotter(x, y, r'$\int_0^{60} WW^v(s)/V_{\ell}(s) ds$' , r'$C(60)/C(0)$', 'cor_C_conc', c_color)
            corr_list.append(corr)
            p_val_list.append(p_val)
        else:
            corr_list.append(0.0)
            p_val_list.append(0.0)
        
        
        
        np.savetxt('corr_list.txt', X=corr_list, delimiter=',',fmt='%1.5f')
        np.savetxt('p_val_list.txt', X=p_val_list, delimiter=',',fmt='%1.5f')
        # scatter plots
        
        
        
        # back to the mother folder
        os.chdir(mother_folder)
        
        
