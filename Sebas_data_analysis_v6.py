# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 21:04:58 2026

@author: Nemat002
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import ast
from collections import deque
import os





class org_class:
    type = 'mixed'
    pass
def matrix_maker(orgs_dict):
    
    w_mat = np.zeros((n_orgs, length))
    c_mat = np.zeros((n_orgs, length))
    
    for key_c in range(len(keys_mixed)):
        key = keys_mixed[key_c]
        
        w_mat[key_c,:] = orgs_dict[key].width_w
        c_mat[key_c,:] = orgs_dict[key].width_c
        
    return w_mat, c_mat

def matrix_plotter(width_w_mat, width_c_mat):

    plt.figure()
    for key_c in range(len(keys_mixed)):
        key = keys_mixed[key_c]
        x = orgs_dict[key].time
        y = width_w_mat[key_c, :]
        plt.plot(x, y, color='m')
        # plt.plot(x, y)
    y_avg = np.mean(width_w_mat, axis=0)
    y_std = np.std(width_w_mat, axis=0)
    y_med = np.median(width_w_mat, axis=0)
    
    plt.plot(x, y_avg, color='k')
    plt.errorbar(x, y_avg, yerr=y_std/np.sqrt(n_orgs-1), capsize=0, color = 'k')
    plt.plot(x, y_med, color='b', label='median')
    plt.xlabel('time (h)')
    plt.ylabel('equivalent '+r'$\ell_W$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('equiv_l_W.PNG', dpi=300)
    
    plt.figure()
    plt.scatter(width_w_mat[:,0], width_w_mat[:,-1])
    
    plt.figure()
    for key_c in range(len(keys_mixed)):
        key = keys_mixed[key_c]
        x = orgs_dict[key].time
        y = width_c_mat[key_c, :]
        plt.plot(x, y, color='g')
        # plt.plot(x, y)
    y_avg = np.mean(width_c_mat, axis=0)
    y_std = np.std(width_c_mat, axis=0)
    y_med = np.median(width_c_mat, axis=0)
    
    plt.plot(x, y_avg, color='k')
    plt.errorbar(x, y_avg, yerr=y_std/np.sqrt(n_orgs-1), capsize=0, color = 'k')
    plt.plot(x, y_med, color='b', label='median')
    plt.xlabel('time (h)')
    plt.ylabel('equivalent '+r'$\ell_C$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('equiv_l_C.PNG', dpi=300)
    
    plt.figure()
    plt.scatter(width_c_mat[:,0], width_c_mat[:,-1])
    
    
    
    return

def plotter(orgs_dict):
    
    plt.figure()
    y_avg = np.zeros(13)
    for key_c in range(len(keys_mixed)):
        key = keys_mixed[key_c]
        x = orgs_dict[key].time
        y = orgs_dict[key].width_w
        plt.plot(x, y, color='m')
        y_avg += y / len(keys_mixed)
    plt.plot(x, y_avg, color='k')
    # plt.plot(x, 0*y_avg+2.0, color='red')
    plt.grid()
    
    # plt.figure()
    # for key_c in range(len(keys_mixed)):
    #     key = keys_mixed[key_c]
    #     x = orgs_dict[key].time
    #     y = orgs_dict[key].width_w
    #     if y[0]>=2:
    #         plt.plot(x, y, color='m')
    
    plt.figure()
    y_avg = np.zeros(13)
    for key_c in range(len(keys_mixed)):
        key = keys_mixed[key_c]
        x = orgs_dict[key].time
        y = orgs_dict[key].width_c
        plt.plot(x, y, color='g')
        y_avg += y / len(keys_mixed)
    plt.plot(x, y_avg, color='k')
    # plt.plot(x, 0*y_avg+2.0, color='red')
    plt.grid()
    
    return

def neigh_mat_maker(neighs_dict, labels, phenos):
    n_cells = len(labels)
    neigh_mat = np.zeros((n_cells, n_cells), dtype=int)
    
    n_neigh            = np.zeros(n_cells, dtype=int)
    n_neigh_diff_pheno = np.zeros(n_cells, dtype=int)
    is_border_sitter   = np.zeros(n_cells, dtype=int)
    
    for i_idx in range(n_cells):
        i_label = labels[i_idx]
        for j_idx in range(i_idx):
            j_label = labels[j_idx]
            
            if j_label in neighs_dict[i_label]:
                # neigh_mat[i_idx, j_idx] = 1
                neigh_mat[i_idx, j_idx] =  ( (-1)**phenos[i_idx] )* ( (-1)**phenos[j_idx] )
                neigh_mat[j_idx, i_idx] =  neigh_mat[i_idx, j_idx]
                # this gives 1 if they are of the same pheno, and -1 if different pheno
    
        
    for i_idx in range(n_cells):    
        n_neigh[i_idx] = int(np.sum(np.abs(neigh_mat[i_idx, :])))
        n_neigh_diff_pheno[i_idx]   = int(np.sum(1*(neigh_mat[i_idx, :]<-0.5)))
        is_border_sitter[i_idx]     = int(bool(n_neigh_diff_pheno[i_idx]))
        
        
    # #sanity check
    # for i in range(n_cells):
    #     indices = list(np.where(neigh_mat[i,:]==1)[0])
    #     neigh_labels = [labels[j] for j in indices]
    #     unity = (set(neighs_dict[labels[i]]) == set(neigh_labels))
    #     print(unity)
        
    return neigh_mat, n_neigh, is_border_sitter , n_neigh_diff_pheno

def shortest_distance_finder(sd_diff_pheno, max_ell_desired , is_border_sitter, neigh_mat):
    
    print("##############################")
    print("Function shortest_distance_finder has a problem! Dont use it. ")
    print("##############################")
    sdfsdfgsdfg
    
    max_ell_current = np.max(sd_diff_pheno)
    if max_ell_current < 0.5: # I mean if it is zero
        sd_diff_pheno[is_border_sitter==1] = 1
        max_ell_current = np.max(sd_diff_pheno)
    
    n_cells = len(sd_diff_pheno)
    dist_temp = max_ell_current
    while dist_temp < max_ell_desired:
        
        for cell_c in range(n_cells):
            
            if cell_c==1:
                dd=4
            if sd_diff_pheno[cell_c]==dist_temp:
                check_list = np.where((neigh_mat[cell_c, :] == 1) & (sd_diff_pheno==0))[0]
                for cell_c_dum in check_list:
                    sd_diff_pheno[cell_c_dum] = dist_temp + 1
        
        dist_temp = dist_temp+1
    
    return sd_diff_pheno


def plot_geometric_graph(x, y, conn, node_types=None, ax=None):
    x = np.asarray(x)
    y = np.asarray(y)
    conn = np.asarray(conn)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    n = len(x)

    # Draw edges
    for i in range(n):
        for j in range(i + 1, n):
            if conn[i, j]:
                ax.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    linewidth=1,
                    alpha=0.5,
                    color="gray"
                )

    # Draw nodes
    if node_types is None:
        ax.scatter(x, y, s=80, zorder=3)
    else:
        node_types = np.asarray(node_types)

        mask1 = node_types == 1
        mask2 = node_types == 0

        ax.scatter(
            x[mask1], y[mask1],
            s=80,
            color="tab:blue",
            label="type 1",
            zorder=3
        )

        ax.scatter(
            x[mask2], y[mask2],
            s=80,
            color="tab:orange",
            label="type 2",
            zorder=3
        )

        # ax.legend()

    # Write node indices
    for i in range(n):
        ax.text(
            x[i],
            y[i],
            str(i),
            fontsize=9,
            ha='center',
            va='center',
            color='k',
            weight='bold',
            zorder=4
        )

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.3)

    return ax

def shortest_distance_matrix(neigh_mat, pheno_diff_mat):
    """
    Computes the shortest path distance between all pairs of nodes
    in an unweighted graph.

    Parameters
    ----------
    connectivity : array-like, shape (n, n)
        0/1 adjacency matrix. connectivity[i][j] = 1 means
        an edge from node i to node j.

    Returns
    -------
    dist : np.ndarray, shape (n, n)
        dist[i][j] is the shortest distance from node i to node j.
        If j is unreachable from i, dist[i][j] = np.inf.
    """
    connectivity = np.abs(neigh_mat)
    A = np.asarray(connectivity)
    n = A.shape[0]

    if A.shape[0] != A.shape[1]:
        raise ValueError("Connectivity matrix must be square.")

    # Convert matrix to adjacency list for efficiency
    neighbors = [np.where(A[i] != 0)[0] for i in range(n)]

    dist = np.full((n, n), np.inf)

    for source in range(n):
        dist[source, source] = 0
        queue = deque([source])

        while queue:
            current = queue.popleft()

            for neighbor in neighbors[current]:
                if np.isinf(dist[source, neighbor]):
                    dist[source, neighbor] = dist[source, current] + 1
                    queue.append(neighbor)
    
    signed_dist = dist * pheno_diff_mat
    
    return signed_dist

def dict_symmetrizer(neighs_dict, labels, label_to_id):
    
    n_cells = len(labels)
    for cell_c_1 in range(n_cells):
        label_1 = int(labels[cell_c_1])
        list_neighs_labels = neighs_dict[label_1]
        
        for label_2 in list_neighs_labels:
            if not (label_1 in neighs_dict[label_2]):
                neighs_dict[label_2].append(label_1)
            
    return neighs_dict

# # testing
# L = 12
# D = 6
# L_0 = 2
# N1 = 8
# N2 = 9
# N_b1 = 2
# N_b2 = 2
# N_tot = N1 + N2


# x = 0.0*np.zeros(N_tot)
# y = 0.0*np.zeros(N_tot)

# types = np.zeros(N_tot, dtype=int)
# types[0:N_b1] = 1
# types[N_b1:N_b1+N_b2] = 0
# types[N_b1+N_b2:(N_b1+N_b2+N1-N_b1)] = 1
# types[(N_b1+N_b2+N1-N_b1):] = 0

# x[0:N_b1] = np.random.uniform(0,L_0, N_b1)
# x[N_b1:N_b1+N_b2] = np.random.uniform(-L_0,0, N_b2)
# x[N_b1+N_b2:(N_b1+N_b2+N1-N_b1)] = np.random.uniform(L_0,L, N1-N_b1)
# x[(N_b1+N_b2+N1-N_b1):] = np.random.uniform(-L, -L_0, N2-N_b2)

# y = np.random.uniform(0,D, N_tot)

# neigh_mat = np.zeros((N_tot, N_tot), dtype=int)
# is_border_sitter = np.zeros(N_tot, dtype=int)
# for i in range(N_tot):
#     for j in range(i):
        
#         r2 = (x[i]-x[j])**2 + (y[i]-y[j])**2
#         if r2 < L_0:
#             neigh_mat[i,j] = 1
#         elif r2 > 15 * L_0:
#             neigh_mat[i,j] = 0
#         else:
#             neigh_mat[i,j] = int(np.random.choice([0,1]))
        
#         if np.abs(r2) < 1e-8:
#             neigh_mat[i,j] = 0
#         neigh_mat[i,j] = neigh_mat[i,j] * (-1)**np.abs(types[i]-types[j])
#         neigh_mat[j,i] = neigh_mat[i,j]

# for i in range(N_tot):    
#     if np.min(neigh_mat[i,:])<-0.5:
#         is_border_sitter[i] = 1
        
# plot_geometric_graph(x, y, np.abs(neigh_mat), types)




# sd = np.zeros(N_tot, dtype=int)
# max_ell_desired = 5
# sd = shortest_distance_finder(sd, max_ell_desired , is_border_sitter, neigh_mat)

# pheno_diff_mat = 1 - 2 * np.bitwise_xor.outer(types, types)
# signed_dist_mat = shortest_distance_matrix(neigh_mat, pheno_diff_mat)

# zdczd
# # testing

knn = int(np.loadtxt('knn.txt', delimiter=','))

if knn==5:
    cols_to_read = ["sample_name", "label", "time", "phenotype_crc_vs_wt", "5_knn_neighbors", "touching_neighbors", "cell_id"]
elif knn==7:
    # cols_to_read = ["sample_name", "label", "time", "phenotype_crc_vs_wt", "7_knn_neighbors", "touching_neighbors", "cell_id"]
    cols_to_read = ["sample_name", "label", "time", "phenotype_crc_vs_wt", "7_knn_neighbors", "touching_neighbors"]
    
# all_data = pd.read_excel('mixed_full_data_test.xlsx', sheet_name='mixed', usecols=cols_to_read)
all_data = pd.read_excel('mixed_70.xlsx', sheet_name='mixed', usecols=cols_to_read)

organoids_names = list(set(list(all_data['sample_name'])))
time_points     = list(set(list(all_data['time'])))
time_points.sort()

n_orgs = len(organoids_names)
n_time_points = len(time_points)

lw_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
# lc_list = [1,2,3,4,5,6,7, 20]
# lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
lc_list = [1,2,3,4,5,6,7, 8, 9, 10, 11, 15, 20, 40]
for l_w in lw_list:
    for l_c in lc_list:

        # l_w = 2
        # l_c = 1
        
        data_folder = 'lw_'+str(int(l_w))+'__lc_'+str(int(l_c))
        if os.path.isdir(data_folder):
            continue
        
        try:
            os.mkdir(data_folder)
        except FileExistsError:
            pass

        conc_switch = 0
        
        shortest_distance_matrices = defaultdict(dict)
        
        # np.full((n_orgs, n_time_points), np.nan)
        w_mat = np.full((n_orgs, n_time_points), np.nan) # population of wt, in each org, at each time; each row in an org;
        c_mat = np.full((n_orgs, n_time_points), np.nan)
        w_a_mat = np.full((n_orgs, n_time_points), np.nan) # affected
        c_a_mat = np.full((n_orgs, n_time_points), np.nan)
        w_v_mat = np.full((n_orgs, n_time_points), np.nan) # visible
        c_v_mat = np.full((n_orgs, n_time_points), np.nan)
        w_b_mat = np.full((n_orgs, n_time_points), np.nan) # border
        c_b_mat = np.full((n_orgs, n_time_points), np.nan)
        
        avg_num_w_visib_to_aff_c_mat = np.full((n_orgs, n_time_points), np.nan) # individual neighberhoods
        avg_num_c_visib_to_aff_w_mat = np.full((n_orgs, n_time_points), np.nan)
        
        
        for org_c in range(n_orgs):
            org_name = organoids_names[org_c]
            print('lw: '+str(l_w)+', '+'lc: '+str(l_c)+', '+'org: '+org_name)
            print('_______________')
            for t_c in range(n_time_points):
                time_val = time_points[t_c]
                
                subset_df = all_data[(all_data["sample_name"] == org_name) & (all_data["time"] == time_val)]
                
                if subset_df.empty: # some organoids do not have 70 h or 65 h data
                    continue
                
                # ids =    np.array(list(subset_df['cell_id']))
                labels = np.array(list(subset_df['label']))
                phenos_str = list(subset_df['phenotype_crc_vs_wt']) # 0 for wt, 1 for cancer
                phenos = np.array([0 if x == 'wt' else 1 for x in phenos_str])
                
                # n_cells = len(ids)
                n_cells = len(labels)
                
                
                # pheno_diff_mat = np.zeros((n_cells, n_cells), dtype=int)
                # for ii in range(n_cells):
                #     for jj in range(n_cells):
                #         if phenos[ii]==phenos[jj]:
                #             pheno_diff_mat[ii,jj] =  1
                #         else:
                #             pheno_diff_mat[ii,jj] = -1
                
                # pheno_diff_mat_2 = (-1)**(np.reshape((1-phenos), (1,len(phenos))).T @ np.reshape(phenos), (1,len(phenos))))
                pheno_diff_mat = 1 - 2 * np.bitwise_xor.outer(phenos, phenos)
                
                
                n_c = np.sum(phenos)
                n_w = n_cells - n_c
                
                # is_sorted = (ids == sorted(ids))
                # if np.any(is_sorted == False):
                #     print("#################")
                #     print("ids is not in sort")
                #     print(org_name)
                #     print(time_val)
                #     print("#################")
                #     sgfdsgfdg
                
                touch_list_list = subset_df["touching_neighbors"].to_list()
                
                if knn==5:
                    knn_list_list   = subset_df["5_knn_neighbors"].to_list()
                elif knn==7:
                    knn_list_list   = subset_df["7_knn_neighbors"].to_list()
                    
                neighs_dict = dict()
                label_to_id = dict()
                for cell_c in range(n_cells):
                    label_to_id[labels[cell_c]] = cell_c
                    if pd.isna(touch_list_list[cell_c]):
                        neighs_candid_touch = []
                    else:
                        neighs_candid_touch = ast.literal_eval(touch_list_list[cell_c])
                    neighs_candid_knn = ast.literal_eval(  knn_list_list[cell_c])
                    
                    if len(neighs_candid_touch)<=len(neighs_candid_knn):    
                        neighs_dict[labels[cell_c]] = neighs_candid_touch.copy()
                    else:
                        neighs_dict[labels[cell_c]] = neighs_candid_knn.copy()
                
                
                neighs_dict = dict_symmetrizer(neighs_dict, labels, label_to_id)
                
                
                neigh_mat, n_neigh, is_border_sitter , n_neigh_diff_pheno = neigh_mat_maker(neighs_dict, labels, phenos)
                # order of rows (and cols) in this matrix is based on the order of 
                # the vector labels (and phenos)
                if not np.all(neigh_mat.T==neigh_mat):
                    print('#######################')
                    print('neigh_mat not symmetric')
                    print('#######################')
                    ffff
                
                n_c_border = np.sum(is_border_sitter * phenos) #exactly on the border
                n_w_border = np.sum(is_border_sitter) - n_c_border #exactly on the border
                
                
                # border_labels_w = # list of their labels
                # border_labels_c = # list of their labels 
                # border_w_cross_neighs = dict()
                # border_c_cross_neighs = dict()
                
                is_affected = np.zeros(n_cells, dtype=int)
                is_visible  = np.zeros(n_cells, dtype=int)
                n_visible_to_self     = np.zeros(n_cells, dtype=int) # how many cells of the other type is visible  to each cell
                n_affected_by_self    = np.zeros(n_cells, dtype=int) # how many cells of the other type is affected by each cell
            
                
                if conc_switch==0:
                    # sd_diff_pheno = np.zeros(n_cells, dtype=int)
                    # shortest distance to the different phenotype
                    signed_dist_mat = shortest_distance_matrix(neigh_mat, pheno_diff_mat)
                    # sd_diff_pheno = shortest_distance_finder(sd_diff_pheno, l_c, is_border_sitter, neigh_mat)
                    
                    # l_c-related
                    cond_lc = (signed_dist_mat < -0.001)     & \
                              (signed_dist_mat > (-l_c-0.001) )
                    lc_related_data = np.where(cond_lc, signed_dist_mat, 0).astype(int)            
                    is_affected_c   = 1*(np.min(lc_related_data, axis=1) * phenos    ).astype(bool)
                    is_visible_w    = 1*(np.min(lc_related_data, axis=1) * (1-phenos)).astype(bool)
                    n_c_affected = np.sum(is_affected_c)
                    n_w_visible  = np.sum(is_visible_w)
                    num_w_neigh_to_c = np.sum(1*(lc_related_data.astype(bool)), axis=1) * (phenos)
                    avg_num_w_visib_to_aff_c =  num_w_neigh_to_c[num_w_neigh_to_c != 0].mean()
                    # l_c-related
                    
                    # l_w-related
                    cond_lw = (signed_dist_mat < -0.001)     & \
                              (signed_dist_mat > (-l_w-0.001) )
                    lw_related_data = np.where(cond_lw, signed_dist_mat, 0).astype(int)
                    is_affected_w   = 1*(np.min(lw_related_data, axis=1) * (1-phenos)).astype(bool)
                    is_visible_c    = 1*(np.min(lw_related_data, axis=1) * phenos    ).astype(bool)
                    n_w_affected = np.sum(is_affected_w)
                    n_c_visible  = np.sum(is_visible_c)
                    num_c_neigh_to_w = np.sum(1*(lw_related_data.astype(bool)), axis=1) * (1-phenos)
                    avg_num_c_visib_to_aff_w =  num_c_neigh_to_w[num_c_neigh_to_w != 0].mean()
                    # l_w-related
                    
                    
                    w_mat[org_c, t_c]   = n_w
                    c_mat[org_c, t_c]   = n_c
                    w_a_mat[org_c, t_c] = n_w_affected
                    c_a_mat[org_c, t_c] = n_c_affected
                    w_v_mat[org_c, t_c] = n_w_visible
                    c_v_mat[org_c, t_c] = n_c_visible
                    w_b_mat[org_c, t_c] = n_w_border
                    c_b_mat[org_c, t_c] = n_c_border
                    
                    avg_num_w_visib_to_aff_c_mat[org_c, t_c] = avg_num_w_visib_to_aff_c
                    avg_num_c_visib_to_aff_w_mat[org_c, t_c] = avg_num_c_visib_to_aff_w
                    
                    
                    
                else:
                    dd
                    # will add this later
                    
                # uyfuyf
                
                # ids.clear()
                # labels.clear()
                # phenos_str.clear()
                # phenos.clear()
                # touch_list_list.clear()
                # knn_list_list.clear()
                # neighs_dict.clear()
        
        
        
        with open(data_folder+"/org_names.txt", "w") as org_names_file:
            for org_name in organoids_names:
                org_names_file.write(org_name + "\n")
        np.savetxt(data_folder+"/time_points.txt", X=time_points, fmt='%.1f', delimiter=',')
        np.savetxt(data_folder+"/lw_lc.txt", X=np.array([l_w, l_c]), fmt='%d', delimiter=',')
        np.savetxt(data_folder+"/w_mat.txt", X=w_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/c_mat.txt", X=c_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/w_a_mat.txt", X=w_a_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/c_a_mat.txt", X=c_a_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/w_v_mat.txt", X=w_v_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/c_v_mat.txt", X=c_v_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/w_b_mat.txt", X=w_b_mat, fmt='%.0f', delimiter=',')
        np.savetxt(data_folder+"/c_b_mat.txt", X=c_b_mat, fmt='%.0f', delimiter=',')
        
        np.savetxt(data_folder+"/avg_num_w_visib_to_aff_c_mat.txt", X=avg_num_w_visib_to_aff_c_mat, fmt='%.3e', delimiter=',')
        np.savetxt(data_folder+"/avg_num_c_visib_to_aff_w_mat.txt", X=avg_num_c_visib_to_aff_w_mat, fmt='%.3e', delimiter=',')






Tamam


n_rows = all_data.shape[0]

keys_mixed = []
orgs_dict = dict()
length = 13
for row_c in range(n_rows):
    
    if (all_data.sample_name[row_c] not in keys_mixed):
        keys_mixed.append(all_data.sample_name[row_c])
        key = all_data.sample_name[row_c]
        orgs_dict[key] = org_class()
        orgs_dict[key].time = np.array(all_data.time[row_c:row_c+length])
        orgs_dict[key].n_w = np.array(all_data.wt_cells[row_c:row_c+length])
        orgs_dict[key].n_c = np.array(all_data.crc_cells[row_c:row_c+length])
        
        orgs_dict[key].n_w_touch = np.array(all_data.interacting_wt_touching[row_c:row_c+length])
        orgs_dict[key].n_c_touch = np.array(all_data.interacting_crc_touching[row_c:row_c+length])
        
        # orgs_dict[key].n_w_touch = np.array(all_data.interacting_wt_13_knn[row_c:row_c+length])
        # orgs_dict[key].n_c_touch = np.array(all_data.interacting_crc_13_knn[row_c:row_c+length])
        
    else:
        continue

n_orgs = len(keys_mixed)
a_w = 100.0
a_c = 121.0
for key_c in range(len(keys_mixed)):
    
    key = keys_mixed[key_c]
    
    orgs_dict[key].A_w = a_w * orgs_dict[key].n_w
    orgs_dict[key].A_c = a_c * orgs_dict[key].n_c
    
    orgs_dict[key].A_w_touch = a_w * orgs_dict[key].n_w_touch
    orgs_dict[key].A_c_touch = a_c * orgs_dict[key].n_c_touch
    
    orgs_dict[key].A_tot = orgs_dict[key].A_w + orgs_dict[key].A_c
    
    orgs_dict[key].r = (orgs_dict[key].A_tot/(4*np.pi))**0.5
    
    orgs_dict[key].th_m = np.arccos( (orgs_dict[key].A_w-orgs_dict[key].A_c) / orgs_dict[key].A_tot)
    
    argument =( (orgs_dict[key].A_w + orgs_dict[key].A_c_touch) \
              - (orgs_dict[key].A_c - orgs_dict[key].A_c_touch) ) \
                / orgs_dict[key].A_tot
                
    orgs_dict[key].th_c = np.arccos( argument )
    
    argument =( (orgs_dict[key].A_w - orgs_dict[key].A_w_touch) \
              - (orgs_dict[key].A_c + orgs_dict[key].A_w_touch) ) \
                / orgs_dict[key].A_tot
                
    orgs_dict[key].th_w = np.arccos( argument )
    
    orgs_dict[key].width_w = ( (orgs_dict[key].th_w - orgs_dict[key].th_m)*orgs_dict[key].r ) / np.sqrt(a_w)
    
    orgs_dict[key].width_c = ( (orgs_dict[key].th_m - orgs_dict[key].th_c)*orgs_dict[key].r ) / np.sqrt(a_c)


# plotter(orgs_dict)
width_w_mat, width_c_mat = matrix_maker(orgs_dict)

matrix_plotter(width_w_mat, width_c_mat)







    
    
    
    
    
    
