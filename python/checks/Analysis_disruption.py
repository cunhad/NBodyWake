#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 11:06:35 2024

@author: asus
"""


# parse stuff

import argparse
import ast



parser = argparse.ArgumentParser(description='disruption anlysis')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--filepath', type=str, help='')
parser.add_argument('--filepath_CP3M', type=str, help='')
parser.add_argument('--sample', type=str, help='')
parser.add_argument('--path_out', type=str, help='')
parser.add_argument('--redshift', type=str, default='3', help='')
parser.add_argument('--redshift_CP3M', type=str, default='3.000', help='')
parser.add_argument('--redshift_CP3M_insertion', type=str, default='10.000', help='')
parser.add_argument('--Nmesh', default='[512,512,32]', help='')
parser.add_argument('--Nmesh_CP3M', default='[2048,2048,2048]', help='')
parser.add_argument('--BoxSize', type=float, default=4, help='')
parser.add_argument('--nfiles', type=int, default=1, help='')
parser.add_argument('--nfiles_CP3M', type=int, default=64, help='')
parser.add_argument('--ncells', type=int, default=2048, help='')
parser.add_argument('--npart', type=int, default=1024, help='')
parser.add_argument('--resol_factor', type=float, default=0.5, help='')
parser.add_argument('--Nside_id', type=int, help='')
parser.add_argument('--do_no_wake', action='store_true', help='')



# parser.add_argument('--num_epochs', type=int, default=10, help='')

args = parser.parse_args()



# parameters

filepath = args.filepath
filepath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/half_lin_cutoff_half_tot_pert_nvpw/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
# filepath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_384/14--11-11pv_1.5708-0.84153-3.0434ra/2dproj/dm/"
# filepath =  '/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_1/-4-11--24pv_0.10211--0.62099-0.7854ra/2dproj/dm/'
# filepath =  '/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE4_tst/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/half_lin_cutoff_half_tot_pert_nvpw/data/1lf_1rf/NSIDE_4/anglid_86/-15--14--17pv_1.4033-0.79983-5.1051ra/2dproj/dm/'
# filepath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/'
# filepath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5029/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_21/-232--108-113pv_0.62237--1.5029-4.4506ra/2dproj/dm/'
# filepath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_86/-153--150--178pv_1.4033-0.79983-5.1051ra/2dproj/dm/'
print("File Path in = "+ str(filepath))

filepath_CP3M = args.filepath_CP3M
filepath_CP3M =  '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/half_lin_cutoff_half_tot_pert_nvpw/'
# filepath_CP3M =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5029/half_lin_cutoff_half_tot_pert_nvpw_v0p6/'
# filepath_CP3M =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/'
print("File Path CUBEP3M in = "+ str(filepath_CP3M))


sample = args.sample
sample =  'sample1001'
# sample =  'sample5001'
# sample =  'sample5029'
print("sample = "+ str(sample))

path_out = args.path_out
path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
# path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE4_tst/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
# path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512/plots_4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4_tst/plots_1/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/plots_3/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
print("Path out = "+ str(path_out))

redshift = args.redshift
redshift = '63'
# redshift = '5'
# redshift = '10'
# redshift = '3'
print("redshift= "+ str(redshift))



redshift_CP3M = args.redshift_CP3M
redshift_CP3M = '63.000'
# redshift_CP3M = '5.000'
# redshift_CP3M = '10.000'
# redshift_CP3M = '3.000'
print("redshift_CP3M= "+ str(redshift_CP3M))



redshift_CP3M_insertion = args.redshift_CP3M_insertion
redshift_CP3M_insertion = '63.000'
# redshift_CP3M_insertion = '5.000'
# redshift_CP3M_insertion = '10.000'
# redshift_CP3M_insertion = '3.000'
print("redshift_CP3M_insertion = "+ str(redshift_CP3M_insertion))





Nmesh =  ast.literal_eval(args.Nmesh)
# Nmesh = [48,48,48]
Nmesh = [48,48,12]
# Nmesh = [512,512,32]
print("Nmesh= "+ str(Nmesh))

Nmesh_CP3M =  ast.literal_eval(args.Nmesh_CP3M)
Nmesh_CP3M = [96,96,96]
print("Nmesh_CP3M= "+ str(Nmesh_CP3M))



BoxSize = args.BoxSize
print("BoxSize= "+ str(BoxSize))

nfiles = args.nfiles
print("nfiles= "+ str(nfiles))

nfiles_CP3M = args.nfiles_CP3M
nfiles_CP3M = 8
print("nfiles_CP3M= "+ str(nfiles_CP3M))

ncells = args.ncells
ncells = 96
print("ncells= "+ str(ncells))

npart = args.npart
npart = 48
print("npart= "+ str(npart))

resol_factor = args.resol_factor
resol_factor = 1
# resol_factor = 2
print("resol_factor= "+ str(resol_factor))

Nside_id = args.Nside_id
# Nside_id = 384
# Nside_id = 1
# Nside_id = 21
Nside_id = 86
print("Nside_id= "+ str(Nside_id))

do_no_wake = args.do_no_wake
do_no_wake = True
print("do_no_wake = "+ str(do_no_wake))

# num_epochs=args.num_epochs
# print("Num epochs = "+ str(num_epochs))





#%%

# import modules


import sys
import os



path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'read')
import Read_slices




#%%

# Read the interpolation files


mesh = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)

if do_no_wake:
    
    
    filepath_nowake = Read_slices.filepath_from_wake_to_nowake(filepath)
    
    mesh_nowake = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath_nowake,redshift)




#%%

# sys.path.append(path_analy+'2d')
# import Projection2d


# Projection2d.plot_2d_proj(mesh)

#%%

# plot the wake particles alongside density contrast

# extract the position of the wake particles

# print("until here2")


sys.path.append(path_analy+'wake_disruption')
import PIDs


pid_wake, _ =  PIDs.extract_pid_wake(Nmesh_CP3M,nfiles_CP3M,ncells,filepath_CP3M,redshift_CP3M_insertion,npart)


# pos wake particles in the CUBEP3M grid size
pos_wake = PIDs.extract_pos_wake(Nmesh_CP3M,nfiles_CP3M,ncells,filepath_CP3M,redshift_CP3M,npart,pid_wake)






#%%


sys.path.append(path_analy+'2d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
import Projection2d



pv, ra = Read_slices.extract_pv_ra(filepath)

# pos wake particles in the Nmesh grid size
depth = Nmesh[2]
pos_wake_rot = PIDs.rotate_pos_wake(pos_wake, ra, pv, ncells, npart, resol_factor,depth)


##nope

# pos_wake_rot[:, [1, 2]] = pos_wake_rot[:, [2, 1]]
# pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [1, 0, 2]]
# pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [2, 0, 1]]


# pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [0, 2, 1]]
# print("until here")


# grid_points_wake = PIDs.obtain_wake_grid_points(Nmesh,pos_wake_rot)

grid_points_wake = PIDs.obtain_wake_grid_points_chunks(Nmesh,pos_wake_rot,chunk_size=pos_wake[:,0].size)

#%%

if Nside_id is not None:
    ns_id = '_nsid' + str(Nside_id)
else:
    ns_id = ''

save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift+ns_id+'wid.png'
PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save_plot_2d_proj_fig)
# # PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save=None)
save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift+ns_id+'.png'
Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)

# print("until here?")

if do_no_wake:
    
    save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift+ns_id+'widd.png'
    PIDs.plot_2d_proj_wakediff_colInfo3d(mesh,mesh_nowake,grid_points_wake,save_plot_2d_proj_fig)





#%%

# plot the wake particles alongside density contrast, for each slice



# grid_points_wake_eachslice = PIDs.obtain_wake_grid_points_eachslice(Nmesh,pos_wake_rot)

grid_points_wake_eachslice = PIDs.obtain_wake_grid_points_eachslice_chunks(Nmesh,pos_wake_rot,chunk_size=pos_wake[:,2].size)


dept = Nmesh[2]
slice_list = list(range(dept))
# slice_list = [1,2]

sliceId =  ["_sl" + str(i)  for i in slice_list]


grid_points_wake_eachslice_list = [grid_points_wake_eachslice[i] for i in slice_list]



# PIDs.plot_2d_proj_wake_colInfo3d_eachSlice(mesh[slice_list], grid_points_wake_eachslice, slice_list,save=None)
save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift+ns_id+sliceId[i]+'wid.png' for i,_ in enumerate(slice_list)]
PIDs.plot_2d_proj_wake_colInfo3d_eachSlice(mesh[:,:,slice_list], grid_points_wake_eachslice, slice_list,save_plot_2d_proj_fig_list)

save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift+ns_id+sliceId[i]+'.png' for i,_ in enumerate(slice_list)]
Projection2d.plot_2d_proj_eachSlice(mesh[:,:,slice_list], slice_list,save_plot_2d_proj_fig_list)


if do_no_wake:
    

    
    save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift+ns_id+sliceId[i]+'widd.png' for i,_ in enumerate(slice_list)]
    PIDs.plot_2d_proj_wakediff_colInfo3d_eachSlice(mesh[:,:,slice_list],mesh_nowake[:,:,slice_list], grid_points_wake_eachslice, slice_list,save_plot_2d_proj_fig_list)



#%%
# sys.path.append(path_analy+'pycurvelab')
# import pycurvelab


slice_list = [1]
sliceId =  ["_sl" + str(i)  for i in slice_list]
grid_points_wake_eachslice_list = [grid_points_wake_eachslice[i] for i in slice_list]
save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift+ns_id+sliceId[i]+'widdc.png' for i,_ in enumerate(slice_list)]
PIDs.plot_2d_proj_curveletFilt_eachSlice(mesh[:,:,slice_list],mesh_nowake[:,:,slice_list], grid_points_wake_eachslice, slice_list,save_plot_2d_proj_fig_list)
# PIDs.plot_2d_proj_curveletFilt_eachSlice(mesh[:,:,slice_list],mesh_nowake[:,:,slice_list], grid_points_wake_eachslice, slice_list)




# pos_wake_rot, shift, lim,Pos = PIDs.rotate_pos_wake(pos_wake, ra, pv, ncells, npart, resol_factor,depth)


#%%



# # import numpy as np

# # n = 48  # Change this to your desired value for n
# # data_tst = np.zeros((n, 3), dtype=np.float32)
# # # Fill the first column with values from 0 to n-1
# # data_tst[:, 0] = np.arange(n, dtype=np.float32)
# # # Fill the third column with 45
# # data_tst[:, 2] = 45.0

# import numpy as np

# n = 48  # Change this to your desired value for n
# m = 48   # Change this to your desired value for m

# # Create an array of shape (n * (m+1), 3) filled with zeros
# data_tst = np.zeros((n * (m+1), 3), dtype=np.float32)

# # Fill the first column with values from 0 to n-1 repeated m+1 times
# data_tst[:, 0] = np.tile(np.arange(n, dtype=np.float32), m+1)

# # Fill the second column with values from 0 to m, repeated n times each
# data_tst[:, 1] = np.repeat(np.arange(m+1, dtype=np.float32), n)

# # Fill the third column with 45
# data_tst[:, 2] = 45.0




#%%

 
 
# Pos = data_tst
# rot_angle = ra
# pivot = pv
# nc = ncells
# nup = npart
# resol_factor = resol_factor
# depth = depth
# lenght_factor=1


# import numpy as np

# # Ensure pivot is a NumPy array to handle element-wise operations
# pivot = np.array(pivot)

# phi, theta, psi  = rot_angle

# Pos = np.mod(Pos, nc)
# Pos = Pos * (nup * resol_factor) / nc

# axis_size = np.array([
#     [nup * resol_factor, 0, 0],
#     [0, nup * resol_factor, 0],
#     [0, 0, nup * resol_factor]
# ])

# # Subtract from Pos, ensuring that pivot is a NumPy array
# Pos -= (nup * resol_factor / 2) + pivot * (nup * resol_factor / nc)

# Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
# Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
# Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

# R = Rz @ Ry @ Rx
# Pos = Pos @ R.T

# axis_size = axis_size @ R.T

# Pos += (1 / (2 * lenght_factor)) * nup * resol_factor

# lim = (1 / lenght_factor) * nup * resol_factor
# Pos_expand = np.empty((0, 3))

# D = np.array([[Dx, Dy, Dz] for Dx in range(-1, 2) for Dy in range(-1, 2) for Dz in range(-1, 2)])
# # D = np.array([[Dx, Dy, Dz] for Dx in range(0, 1) for Dy in range(0, 1) for Dz in range(0, 1)])
# shifts = D @ axis_size

# for shift in shifts:
#     Pos_aux = Pos + shift
#     mask = np.all((Pos_aux >= 0) & (Pos_aux < lim), axis=1)
#     Pos_expand = np.vstack([Pos_expand, Pos_aux[mask]])
#     # mask = np.all((Pos_aux >= 0) & (Pos_aux <= lim))
#     # Pos_expand = np.vstack([Pos_expand, Pos_aux[mask]])
#     # Pos_expand = np.vstack([Pos_expand, Pos_aux])
    
# Pos_expand[:,0:2] -= 0.5  
# Pos_expand[:,2] = Pos_expand[:,2]/(lim/depth)

# Pos[:,2] = Pos[:,2]/(lim/depth)

# # return Pos_expand, shift, lim, Pos

# #%%

# # grid_points_data_tst = PIDs.obtain_wake_grid_points_chunks(Nmesh,data_tst,chunk_size=pos_wake[:,0].size)
# # grid_points_data_tst_eachslice = PIDs.obtain_wake_grid_points_eachslice_chunks(Nmesh,data_tst,chunk_size=pos_wake[:,2].size)


# data_tst_rot = PIDs.rotate_pos_wake(data_tst, ra, pv, ncells, npart, resol_factor,depth)
# grid_points_data_tst_rot = PIDs.obtain_wake_grid_points_chunks(Nmesh,data_tst_rot,chunk_size=pos_wake[:,0].size)
# grid_points_data_tst_rot_eachslice = PIDs.obtain_wake_grid_points_eachslice_chunks(Nmesh,grid_points_data_tst_rot,chunk_size=pos_wake[:,2].size)



# #%%

# PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_data_tst_rot,save_plot_2d_proj_fig)


# PIDs.plot_2d_proj_wake_colInfo3d_eachSlice(mesh[:,:,slice_list], grid_points_data_tst_rot_eachslice, slice_list,save_plot_2d_proj_fig_list)



#%%

# #%%

# import numpy as np

# Pos =  np.array([[24,1,47],[1,1,48]])
# rot_angle = ra
# pivot = pv
# nc = ncells
# nup = npart
# resol_factor = 1
# lenght_factor=1






# # Ensure pivot is a NumPy array to handle element-wise operations
# pivot = np.array(pivot)

# phi, theta, psi = rot_angle

# Pos = np.mod(Pos, nc)
# Pos = Pos * (nup * resol_factor) / nc

# axis_size = np.array([
#     [nup * resol_factor, 0, 0],
#     [0, nup * resol_factor, 0],
#     [0, 0, nup * resol_factor]
# ])

# # Subtract from Pos, ensuring that pivot is a NumPy array
# Pos -= (nup * resol_factor / 2) + pivot * (nup * resol_factor / nc)

# Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
# Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
# Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

# R = Rz @ Ry @ Rx
# Pos = Pos @ R.T

# axis_size = axis_size @ R.T

# Pos += (1 / (2 * lenght_factor)) * nup * resol_factor

# lim = (1 / lenght_factor) * nup * resol_factor
# Pos_expand = np.empty((0, 3))

# D = np.array([[Dx, Dy, Dz] for Dx in range(-1, 2) for Dy in range(-1, 2) for Dz in range(-1, 2)])
# shifts = D @ axis_size

# for shift in shifts:
#     Pos_aux = Pos + shift
#     mask = np.all((Pos_aux >= 0) & (Pos_aux <= lim), axis=1)
#     Pos_expand = np.vstack([Pos_expand, Pos_aux[mask]])



