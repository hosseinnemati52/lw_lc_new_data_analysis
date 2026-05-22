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

lw_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]

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

for l_w in lw_list:
    for l_c in lc_list:
        
        # go to daughter folder
        folder_name = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        try:
            os.chdir(mother_folder + "\\"+ folder_name) #windows
        except:
            os.chdir(mother_folder + "/"+ folder_name) #linux
        
        time_points = np.loadtxt("time_points.txt", delimiter=',')
        
        w_mat = np.loadtxt("w_mat.txt", delimiter=',')
        w_b_mat = np.loadtxt("w_b_mat.txt", delimiter=',')
        w_a_mat = np.loadtxt("w_a_mat.txt", delimiter=',')
        w_v_mat = np.loadtxt("w_v_mat.txt", delimiter=',')
        w_u_mat = w_mat - w_a_mat
        
        c_mat = np.loadtxt("c_mat.txt", delimiter=',')
        c_b_mat = np.loadtxt("c_b_mat.txt", delimiter=',')
        c_a_mat = np.loadtxt("c_a_mat.txt", delimiter=',')
        c_v_mat = np.loadtxt("c_v_mat.txt", delimiter=',')
        c_u_mat = c_mat - c_a_mat
        
        np.savetxt("c_u_mat.txt", X=c_u_mat, fmt='%d', delimiter=',')
        np.savetxt("w_u_mat.txt", X=w_u_mat, fmt='%d', delimiter=',')
        
        dt = time_points[1] - time_points[0]
        
        n_orgs = np.shape(w_mat)[0]
        n_time = np.shape(w_mat)[1]
        
        beta_w_a = np.zeros((n_orgs, n_time))
        beta_c_a = np.zeros((n_orgs, n_time))
        
        beta_w_a[:,0] = beta_w_0*0.5
        beta_c_a[:,0] = beta_c_0*0.5
        
        vol_lum = np.zeros((n_orgs, n_time))
        
        fw_mat = w_a_mat/w_mat # fraction of affected
        fc_mat = c_a_mat/c_mat
        
        fw_beta_a_w = np.zeros((n_orgs, n_time))
        fc_beta_a_c = np.zeros((n_orgs, n_time))
        
        fw_beta_a_w[:,0] = beta_w_a[:,0] * fw_mat[:,0]
        fc_beta_a_c[:,0] = beta_c_a[:,0] * fc_mat[:,0]
        
        fw_beta_a_w_integ = np.zeros((n_orgs, n_time))
        fc_beta_a_c_integ = np.zeros((n_orgs, n_time))
        
        Wa_beta_a_w_integ = np.zeros((n_orgs, n_time))
        Ca_beta_a_c_integ = np.zeros((n_orgs, n_time))

        # dlogWdt = np.zeros((n_orgs, n_time))
        # dlogCdt = np.zeros((n_orgs, n_time))
        
        logWtoW0 = np.log(w_mat/w_mat[:,[0]])
        logCtoC0 = np.log(c_mat/c_mat[:,[0]])
        
        # dlogWdt[:, 0] = (np.log(w_mat[:, 1])-np.log(w_mat[:, 0]))/dt
        # dlogCdt[:, 0] = (np.log(c_mat[:, 1])-np.log(c_mat[:, 0]))/dt
        
        # dlogWdt[:, -1] = (np.log(w_mat[:, -1])-np.log(w_mat[:, -2]))/dt
        # dlogCdt[:, -1] = (np.log(c_mat[:, -1])-np.log(c_mat[:, -2]))/dt
        
        for org_c in range(n_orgs):
            for t_c in range(0, len(time_points)):
                
                v_shell = v_c_avg_w * w_mat[org_c, t_c] + v_c_avg_c * c_mat[org_c, t_c]
                A = 3 * shell_thickness 
                B = 3 * shell_thickness**2
                C = shell_thickness**3 - 3*v_shell/(4 * np.pi)
                delta = B**2 - 4 * A * C
                if delta <0:
                    print('negative delta')
                    shfshf
                R_i = (-B + np.sqrt(delta))/(2*A)
                vol_lum[org_c, t_c] = (4/3)*np.pi*(R_i**3)
                
                # if t_c>0 and t_c<len(time_points)-1:
                if t_c>0:
                    # beta_w_a[org_c, t_c] = (1/w_a_mat[org_c, t_c]) * \
                    #                           ( (2/dt)*(w_mat[org_c, t_c]-w_mat[org_c, t_c-1])  \
                    #                            - beta_w_0 * (w_u_mat[org_c, t_c-1]+w_u_mat[org_c, t_c]) \
                    #                            - beta_w_a[org_c, t_c-1] * w_a_mat[org_c, t_c-1] )
                                                  
                    # beta_c_a[org_c, t_c] = (1/c_a_mat[org_c, t_c]) * \
                    #                           ( (2/dt)*(c_mat[org_c, t_c]-c_mat[org_c, t_c-1])  \
                    #                            - beta_c_0 * (c_u_mat[org_c, t_c-1]+c_u_mat[org_c, t_c]) \
                    #                            - beta_c_a[org_c, t_c-1] * c_a_mat[org_c, t_c-1] )
                    
                    # dlogWdt[org_c, t_c] = (np.log(w_mat[org_c, t_c+1])-np.log(w_mat[org_c, t_c-1]))/(2*dt)
                    # dlogCdt[org_c, t_c] = (np.log(c_mat[org_c, t_c+1])-np.log(c_mat[org_c, t_c-1]))/(2*dt)
                    
           
                    fw_integ = np.sum(fw_mat[org_c,1:t_c])*dt + 0.5*dt*(fw_mat[org_c,0]+fw_mat[org_c,t_c])
                    fc_integ = np.sum(fc_mat[org_c,1:t_c])*dt + 0.5*dt*(fc_mat[org_c,0]+fc_mat[org_c,t_c])
                    
                    fw_beta_a_w_integ[org_c, t_c] = logWtoW0[org_c, t_c] - beta_w_0*time_points[t_c] + beta_w_0 * fw_integ
                    fc_beta_a_c_integ[org_c, t_c] = logCtoC0[org_c, t_c] - beta_c_0*time_points[t_c] + beta_c_0 * fc_integ
                    
                    
                    w_u_integ = np.sum(w_u_mat[org_c,1:t_c])*dt + 0.5*dt*(w_u_mat[org_c,0]+w_u_mat[org_c,t_c])
                    c_u_integ = np.sum(c_u_mat[org_c,1:t_c])*dt + 0.5*dt*(c_u_mat[org_c,0]+c_u_mat[org_c,t_c])
                    
                    Wa_beta_a_w_integ[org_c, t_c] = (w_mat[org_c, t_c]-w_mat[org_c, 0]) - beta_w_0*w_u_integ
                    Ca_beta_a_c_integ[org_c, t_c] = (c_mat[org_c, t_c]-c_mat[org_c, 0]) - beta_c_0*c_u_integ
                    
                    # fw_beta_a_w[org_c, t_c] = (fw_beta_a_w_integ[org_c, t_c] - fw_beta_a_w_integ[org_c, t_c-1]) *(2/dt) - fw_beta_a_w[org_c, t_c-1]
                    # fc_beta_a_c[org_c, t_c] = (fc_beta_a_c_integ[org_c, t_c] - fc_beta_a_c_integ[org_c, t_c-1]) *(2/dt) - fc_beta_a_c[org_c, t_c-1]
                    
                # beta_w_a[org_c, t_c] = dlogWdt[org_c, t_c] / fw_mat[org_c, t_c] + beta_w_0 * (1-1/fw_mat[org_c, t_c])
                # beta_c_a[org_c, t_c] = dlogCdt[org_c, t_c] / fc_mat[org_c, t_c] + beta_c_0 * (1-1/fc_mat[org_c, t_c])
                
                
        # fsf
        
        # np.savetxt("beta_w_a.txt", X=beta_w_a, fmt='%.4e', delimiter=',')
        # np.savetxt("beta_c_a.txt", X=beta_w_a, fmt='%.4e', delimiter=',')
        
        np.savetxt("vol_lum.txt", X=vol_lum, fmt='%.4e', delimiter=',')
        
        
        # np.savetxt("fw_beta_a_w.txt", X=fw_beta_a_w, fmt='%.4e', delimiter=',')
        # np.savetxt("fc_beta_a_c.txt", X=fc_beta_a_c, fmt='%.4e', delimiter=',')
        
        np.savetxt("fw_beta_a_w_integ.txt", X=fw_beta_a_w_integ, fmt='%.4e', delimiter=',')
        np.savetxt("fc_beta_a_c_integ.txt", X=fc_beta_a_c_integ, fmt='%.4e', delimiter=',')
        
        np.savetxt("Wa_beta_a_w_integ.txt", X=Wa_beta_a_w_integ, fmt='%.4e', delimiter=',')
        np.savetxt("Ca_beta_a_c_integ.txt", X=Ca_beta_a_c_integ, fmt='%.4e', delimiter=',')
        
        
        # plt.figure()
        # for org_c in range(n_orgs):
        #     plt.plot(time_points, beta_w_a[org_c, :], color='m')
        # plt.plot(time_points, beta_w_0+0.0*time_points , color='k', linestyle='--')
        # plt.plot(time_points, np.mean(beta_w_a, axis=0) , color='r', linestyle='--')
        # plt.ylim(0,1.2*beta_w_0)
        # plt.savefig('beta_w_a.PNG', dpi=300)
        
        # plt.figure()
        # for org_c in range(n_orgs):
        #     plt.plot(time_points, beta_c_a[org_c, :], color='g')
        # plt.plot(time_points, beta_c_0+0.0*time_points , color='k', linestyle='--')
        # plt.plot(time_points, np.mean(beta_c_a, axis=0) , color='r', linestyle='--')
        # plt.ylim(0,1.2*beta_c_0)
        # plt.savefig('beta_c_a.PNG', dpi=300)
        
        
        
        # back to the mother folder
        os.chdir(mother_folder)
        
        