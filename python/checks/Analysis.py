#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:07:14 2023

@author: Disrael
"""


# Input form CUBEP3M
import sys
import os

path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'ToMesh')
import ToMeshCUBEP3M

# filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/'
# filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps48_48/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
redshift = '5.000'
redshift_ = '5'
Nmesh = [12,48,48]
BoxSize = 96
nfiles = 8
ncells = 96
npart = 48
sample = "sample1001"


# # mesh = ToMeshCUBEP3M.readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
mesh, _, _ = ToMeshCUBEP3M.readCUBEP3M2(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
# # mesh = ToMeshCUBEP3M.Mesh_Wake(Nmesh,BoxSize)


#%%
# Input from grid3d Binary files

import sys
import os



path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'read')
import Read_slices


# From grid3d grid3d Binary files at an angle from folder



# path_input =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512/4Mpc_2048c_1024p_zi63_nowakem/sample5001/data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
# path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512/plots/4Mpc_2048c_1024p_zi63_nowakem/"
# sample = "sample5001"

# filepath = path_input
# redshift = '3'
# redshift_ = '3'
# Nmesh = 512
# BoxSize = 4
# nfiles = 32
# depth = 32


path_input =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
sample = "sample1001"

filepath = path_input
redshift = '5'
redshift_ = '5'
# Nmesh = 48
Nmesh = [12,48,48]
BoxSize = 4
nfiles = 1

mesh2 = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)

# pos_wake = PIDs.extract_pos_wake(Nmesh_o,nfiles,ncells,filepath,redshift,npart)


#%%

# # From grid3d grid3d Binary files at an angle anglid helpix



# path_input =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE8/4Mpc_2048c_1024p_zi63_nowakem/sample5046/data/1lf_0.5rf/NSIDE_8/anglid_1/-43-113--256pv_0.10211--0.62099-0.7854ra/2dproj/dm/"
# path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE8/plots2/4Mpc_2048c_1024p_zi63_nowakem/"
# sample = "sample5046"
# anglid = "_aid1"

# filepath = path_input
# redshift = '3'
# redshift_ = '3'
# Nmesh = 512
# BoxSize = 4
# nfiles = 1
# depth = 32

# path_input =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps512_512/4Mpc_2048c_1024p_zi63_nowakem/sample5001/data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
# path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps512_512/plots2/4Mpc_2048c_1024p_zi63_nowakem/"
# sample = "sample5001"
# anglid = '_aid0'

# filepath = path_input
# redshift = '3'
# redshift_ = '3'
# Nmesh = 512
# BoxSize = 4
# nfiles = 1
# depth = 512

# mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles, depth,filepath,redshift)





# path_input =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_1/-4-11--24pv_0.10211--0.62099-0.7854ra/2dproj/dm/"
path_input =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_384/14--11-11pv_1.5708-0.84153-3.0434ra/2dproj/dm/"
sample = "sample1001"

filepath = path_input
redshift = '5'
redshift_ = '5'
# Nmesh = 48
Nmesh = [12,48,48]
BoxSize = 4
nfiles = 1

mesh2 = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)




#%%


# plot the wake particles alongside density contrast

import sys
import os



path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'read')
import Read_slices


# path_input_bin =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
# path_input_bin =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_1/-4-11--24pv_0.10211--0.62099-0.7854ra/2dproj/dm/"
path_input_bin =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_384/14--11-11pv_1.5708-0.84153-3.0434ra/2dproj/dm/"
# path_input =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_1/-4-11--24pv_0.10211--0.62099-0.7854ra/2dproj/dm/"

# extract the pivot and the rotation angles
pv, ra = Read_slices.extract_pv_ra(path_input_bin)

# extract the position of the wake particles

sys.path.append(path_analy+'wake_disruption')
import PIDs

path_input = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/'
filepath = path_input
redshift = '5.000'
redshift_ = '5'
Nmesh_o = [96,96,96]
BoxSize = 4
nfiles = 8
ncells = 96
npart = 48
sample = "sample1001"

pos_wake = PIDs.extract_pos_wake(Nmesh_o,nfiles,ncells,filepath,redshift,npart)


nc = ncells
nup = npart
resol_factor = 1
Pos = pos_wake
pivot = pv
rot_angle = ra


pos_wake_rot = PIDs.rotate_pos_wake(pos_wake, ra, pv, ncells, npart, resol_factor)

pos_wake_rot[:, [1, 2]] = pos_wake_rot[:, [2, 1]]

pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [1, 0, 2]]


#%%

grid_points_wake = PIDs.obtain_wake_grid_points(Nmesh,pos_wake_rot)

PIDs.plot_2d_proj_wake_colInfo3d(mesh2,grid_points_wake,save=None)



#%%

# # Plot Projection

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
sys.path.append(path_analy+'2d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
import Projection2d

# # save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+anglid+'.png'
# # # save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+'.png'
# Projection2d.plot_2d_proj(mesh)
# # proj_2d= Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)

Projection2d.plot_2d_proj(mesh2)


# # # Plot Projection slices

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# sys.path.append(path_analy+'2d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# import Projection2d

# # slice_list = [0]
# # # slice_list = list(range(32))
# # sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]
# # deept = 1

# dept = 4
# slice_list = list(range(int(32/dept)))
# sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]

# save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]
# test= Projection2d.plot_2d_proj_eachSlice(mesh,slice_list,dept,save_plot_2d_proj_fig_list)











# # histogram 3D

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append(path_analy+'3d')

# import Analysis3d

# save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+anglid+'.png'
# Analysis3d.histogram_3d(mesh,save_3Dhist_filename)






# # # histogram 3D slices

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append(path_analy+'3d')
# import Analysis3d

# # # slice_list = [0]
# # slice_list = list(range(32))
# # sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]
# # dept = 1

# dept = 4
# slice_list = list(range(int(32/dept)))
# sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]

# save_3Dhist_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]
# Analysis3d.histogram_3d_eachSlice(mesh,slice_list,dept,save_3Dhist_2d_proj_fig_list)








# # Power Spectrum

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# sys.path.append(path_analy+'ps')
# import powerSpectrum_nbodykit

# save_PS_filename = path_out+'PS_'+sample+'_z'+redshift_+anglid+'.png'
# powerSpectrum_nbodykit.powerSpectrum(mesh,save_PS_filename)




# # compute the 2D power

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# sys.path.append(path_analy+'ps')
# import powerSpectrum_nbodykit

# # nmu=5
# # # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# nmu=5
# # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# save_PS2D_filename = path_out+'PS2D_'+sample+'_z'+redshift_+'.png'
# # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu,save_PS2D_filename)


#%%

# Obtain the PIDs






sys.path.append(path_analy+'wake_disruption')
import PIDs



pos_wake = PIDs.extract_pos_wake(Nmesh,nfiles,ncells,filepath,redshift,npart)

grid_points_wake = PIDs.obtain_wake_grid_points(Nmesh,pos_wake)

PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save=None)




# PIDs.plot_2d_proj_wake(mesh,pos_wake,save=None)







# grid_points_wake_eachSlice = PIDs.obtain_wake_grid_points_eachslice(Nmesh,pos_wake)


# i=11
# PIDs.plot_2d_slice_wake_colInfo2d(mesh[i,:,:],grid_points_wake_eachSlice[i],save=None)










#%%

# grid_points_wake2 = PIDs.obtain_wake_grid_points(Nmesh,pos_wake)


# grid_points_wake3 = PIDs.obtain_wake_grid_points_eachslice(Nmesh,pos_wake)
# PIDs.plot_2d_proj_wake_colInfo(mesh,pos_wake,save=None)
# PIDs.plot_2d_proj_wake_colInfo3d(mesh,pos_wake,save=None)

# grid_points_wake_eachSlice = PIDs.obtain_wake_grid_points_eachslice(Nmesh,pos_wake)

#%%

# PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save=None)


# i=11
# PIDs.plot_2d_slice_wake_colInfo2d(mesh[i,:,:],grid_points_wake_eachSlice[i],save=None)



#%%

# tst = mesh[:,:,0]


# import numpy as np


# # # find the ranges in which the particles in the wake are
# # find the cubes in which the particles of the wake are

# # nfiles = 8
# # npart = 4

# cubes_per_dim = int(round(nfiles ** (1/3)))
# cubes_per_plane = cubes_per_dim * cubes_per_dim
# cubes_below_wake = list(range(int(nfiles / 2 - cubes_per_plane), int(nfiles / 2)))
# cubes_above_wake = list(range(int(nfiles / 2), int(nfiles / 2 + cubes_per_plane)))

# part_wake_id_start = []
# part_wake_id_end = []

# npart_dim_cube = int(npart / cubes_per_dim)
# nun_part_slabcube = int(npart_dim_cube * npart_dim_cube)
# nun_part_cube  = int(npart_dim_cube * npart_dim_cube * npart_dim_cube)

# for cub in cubes_below_wake:
#     part_wake_id_start.append(((cub + 1) * nun_part_cube) - nun_part_slabcube + 1)
#     part_wake_id_end.append((cub + 1) * nun_part_cube)

# for cub in cubes_above_wake:
#     part_wake_id_start.append((cub * nun_part_cube)  + 1)
#     part_wake_id_end.append((cub * nun_part_cube) + nun_part_slabcube)

# # 

# # pid_wake_tot = []
# # xv_wake_tot = []

# # Initialize empty arrays to store combined results
# # pid_wake_tot = np.empty((0, 1))  # Start with an empty array of shape (0, 1)
# pid_wake_tot = np.array([])     # Start with an empty vector (1D array)
# xv_wake_tot = np.empty((0, 3))  # Start with an empty array of shape (0, 3)


# # nod = 1
# # for i in range(nod,nod+1):
# for i in range(0,nfiles):    
    
#     # filename = filepath+redshift+'PID'+str(i)+'.dat'        
#     filename = filepath+'PID'+str(i)+'.ic'        
#     # data_pid = np.fromfile(filename, dtype=np.integer, offset=12*4)
#     data_pid = np.fromfile(filename, dtype=np.integer, offset=4)    #for ic

   
#     # filename = filepath+redshift+'xv'+str(i)+'.dat'
#     filename = filepath+'xv'+str(i)+'.ic'
#     node = int(i)
#     nc = ncells
#     number_node_dim = nfiles**(1./3)   
#     k_node = np.floor(node/number_node_dim**2)
#     res = np.floor(node % number_node_dim**2)
#     j_node = np.floor(res/number_node_dim);
#     i_node=res % number_node_dim
    
#     # data_xv = np.fromfile(filename, dtype=np.float32 , offset=4*(12)).reshape((-1,6))
#     data_xv = np.fromfile(filename, dtype=np.float32 , offset=4).reshape((-1,6)) # for ic
    
#     data_xv[:,0] = data_xv[:,0] + (nc/number_node_dim)*i_node
#     data_xv[:,1] = data_xv[:,1] + (nc/number_node_dim)*j_node
#     data_xv[:,2] = data_xv[:,2] + (nc/number_node_dim)*k_node
    
#     # range_wake = range(np * np * (np - 1) / 2, np * np * (np + 1) / 2)
    
#     # List to store index of elements that are within any of the ranges
#     within_range_elements = []
    
#     # Iterate over each element in lst2
#     for index, value in enumerate(data_pid):
#         # Check if the element is within any range
#         for i in range(len(part_wake_id_start)):
#             if part_wake_id_start[i] <= value <= part_wake_id_end[i]:
#                 within_range_elements.append(index)
#                 break  # Break once we find a range, no need to check further
    
#     # indexes_in_range = [index for index, value in enumerate(data_pid) if (npart * npart * (npart - 1) / 2) <= value <= (npart * npart * (npart + 1) / 2)]

#     pid_wake = data_pid[within_range_elements]
#     xv_wake = data_xv[within_range_elements,0:3]
    
#     # pid_wake_tot.extend(pid_wake)
#     # xv_wake_tot.extend(* xv_wake)
    
#     # Combine with the previous arrays
#     pid_wake_tot = np.concatenate([pid_wake_tot, pid_wake])  # Combines n x 1 arrays into ? x 1
#     xv_wake_tot = np.vstack([xv_wake_tot, xv_wake])  # Combines n x 3 arrays into ? x 3


#     # xv_wake_tot = [combined_list_3xN_row + new_row for combined_list_3xN_row, new_row in zip(xv_wake_tot, xv_wake)]












    
# # data_sorted = sorted(data)
# indexed_numbers = list(enumerate(data_pid))
# sorted_indexed_numbers = sorted(indexed_numbers, key=lambda x: x[1])
# indices, sorted_pid  = zip(*sorted_indexed_numbers)
# sorted_xv = data_xv[indices,:]

# sorted_pid_lst = np.int32(list(sorted_pid))






# nc = ncells
# number_node_dim = nfiles**(1./3)
# nc_dim = 

# selected_indices = [i for i, x in enumerate(list1) if x > 0]
# selected_list1 = [list1[i] for i in selected_indices]
# selected_list2 = [list2[i] for i in selected_indices]
    
