#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 19:03:45 2023

@author: asus
"""

# command line arguments
import sys
import os

path_analy = os.getcwd() +'/' 
# path_input =  sys.argv[1]
# path_out = sys.argv[2]
# sample = sys.argv[3]
# anglid = sys.argv[4]
depth = 32



sys.path.append(path_analy+'read')
import Read_slices

path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_test/"

filepath =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001//half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/"
redshift = '3'
redshift_ = '3'
Nmesh = [512,512,32]
BoxSize = 4
nfiles = 1
# nfiles = 512
sample = "sample5001"
anglid = "_aid1"

# mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)
mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)
#%%

slic = 1

mesh2d= Read_slices.read_slice_bin(Nmesh,BoxSize,filepath,redshift,slic)

#%%
sys.path.append(path_analy+'2d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
import Projection2d


slice_list = [slic]
# slice_list = list(range(32))
sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]

# sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]


save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]


Projection2d.plot_2d_proj_eachSlice(mesh,slice_list,save_plot_2d_proj_fig_list)


#%%

# import matplotlib
# from matplotlib import pyplot as plt
# import numpy as np

# matplotlib.use('Qt5Agg')
# # plt.switch_backend('Qt5Agg')
# # matplotlib.use("TkAgg")  # or "Qt5Agg" if you prefer


# plt.figure()    
# plt.imshow(np.log10(+1+mesh[:,:,0]))

# # plt.imshow(mesh[:,:,0], cmap='viridis')  # you can use other colormaps like 'gray', 'plasma', etc.

#%%


save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'_t.png' for i in slice_list]
Projection2d.plot_2d_proj_eachSlice(mesh2d,slice_list,save_plot_2d_proj_fig_list)



#%%

filename_ = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001//half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin"
import numpy as np

dt = np.fromfile(filename_, dtype=np.float32)






