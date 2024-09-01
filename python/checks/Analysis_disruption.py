#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 11:06:35 2024

@author: asus
"""


# parse stuff

import argparse

parser = argparse.ArgumentParser(description='disruption anlysis')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--filepath', type=str, help='')
parser.add_argument('--filepath_CP3M', type=str, help='')
parser.add_argument('--sample', type=str, help='')
parser.add_argument('--path_out', type=str, help='')
parser.add_argument('--redshift', type=str, default=3, help='')
parser.add_argument('--redshift_CP3M', type=str, default=3.000, help='')
parser.add_argument('--Nmesh', type=list, default=[32,512,512], help='')
parser.add_argument('--Nmesh_CP3M', type=list, default=[512,512,512], help='')
parser.add_argument('--BoxSize', type=float, default=4, help='')
parser.add_argument('--nfiles', type=int, default=1, help='')
parser.add_argument('--nfiles_CP3M', type=int, default=64, help='')
parser.add_argument('--ncells', type=int, default=2048, help='')
parser.add_argument('--npart', type=int, default=1024, help='')
parser.add_argument('--resol_factor', type=float, default=0.5, help='')
parser.add_argument('--Nside_id', type=int, help='')


# parser.add_argument('--num_epochs', type=int, default=10, help='')

args = parser.parse_args()



# parameters

# filepath = args.filepath
# filepath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
filepath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_384/14--11-11pv_1.5708-0.84153-3.0434ra/2dproj/dm/"
print("File Path in = "+ str(filepath))

# filepath_CP3M = args.filepath_CP3M
filepath_CP3M =  '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/'
print("File Path CUBEP3M in = "+ str(filepath_CP3M))


sample = args.sample
sample =  'sample1001'
print("sample = "+ str(sample))

path_out = args.path_out
# path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
path_out =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
print("Path out = "+ str(path_out))

redshift = args.redshift
redshift = '5'
print("Redshift= "+ str(redshift))

redshift_CP3M = args.redshift_CP3M
redshift_CP3M = '5.000'
print("redshift_CP3M= "+ str(redshift_CP3M))


Nmesh = args.Nmesh
Nmesh = [12,48,48]
print("Nmesh= "+ str(Nmesh))

Nmesh_CP3M = args.Nmesh_CP3M
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
print("resol_factor= "+ str(resol_factor))

Nside_id = args.Nside_id
Nside_id = 384
print("Nside_id= "+ str(Nside_id))

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


#%%

# plot the wake particles alongside density contrast

# extract the position of the wake particles

sys.path.append(path_analy+'wake_disruption')
import PIDs


pos_wake = PIDs.extract_pos_wake(Nmesh_CP3M,nfiles_CP3M,ncells,filepath_CP3M,redshift_CP3M,npart)

pv, ra = Read_slices.extract_pv_ra(filepath)

pos_wake_rot = PIDs.rotate_pos_wake(pos_wake, ra, pv, ncells, npart, resol_factor)
# pos_wake_rot[:, [1, 2]] = pos_wake_rot[:, [2, 1]]
# pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [1, 0, 2]]
pos_wake_rot[:, [0, 1, 2]] = pos_wake_rot[:, [2, 0, 1]]

grid_points_wake = PIDs.obtain_wake_grid_points(Nmesh,pos_wake_rot)

if Nside_id is not None:
    ns_id = '_nsid' + str(Nside_id)
else:
    ns_id = ''

# save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift+ns_id+'.png'
# PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save_plot_2d_proj_fig)
# # PIDs.plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save=None)


#%%




grid_points_wake_eachslice = PIDs.obtain_wake_grid_points_eachslice(Nmesh,pos_wake_rot)

dept = Nmesh[0]
slice_list = list(range(dept))
# slice_list = [1,2]

sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]


grid_points_wake_eachslice_list = [grid_points_wake_eachslice[i] for i in slice_list]



# PIDs.plot_2d_proj_wake_colInfo3d_eachSlice(mesh[slice_list], grid_points_wake_eachslice, slice_list,save=None)
save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift+ns_id+sliceId[i]+'.png' for i in slice_list]
PIDs.plot_2d_proj_wake_colInfo3d_eachSlice(mesh[slice_list], grid_points_wake_eachslice, slice_list,save_plot_2d_proj_fig_list)








