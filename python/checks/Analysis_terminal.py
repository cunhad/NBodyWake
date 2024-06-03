#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:07:14 2023

@author: Disrael
"""



# command line arguments
import sys
import os


# # total arguments
# n = len(sys.argv)
# print("Total arguments passed:", n)

# # Arguments passed
# print("\nName of Python script:", sys.argv[0])


# print("\nArguments passed:", end = " ")
# for i in range(1, n):
#     print(sys.argv[i], end = " ")

# example of argument:
    # Input form CUBEP3M:
    # /home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/ /home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/ sample5001 _aid1 32
    # Input from Binary files:
    # /home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512/4Mpc_2048c_1024p_zi63_nowakem/sample5001/data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/   /home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512/plots/4Mpc_2048c_1024p_zi63_nowakem/ sample5001 aid0 32

# print(sys.argv[1].split('/')[-2])
# splited = sys.argv[0].split('/')   
# path_analy = "/".join(splited[0:-1])+'/'
# # path_analy =  sys.argv[1]
path_analy = os.getcwd() +'/' 
path_input =  sys.argv[1]
path_out = sys.argv[2]
sample = sys.argv[3]
anglid = sys.argv[4]
depth = int(sys.argv[5])


# if sys.argv[1].split('/')[-2][0] == 's':
#     sample = sys.argv[1].split('/')[-2]
# else:
#     sample = sys.argv[1].split('/')[-3]
    

# # Input form CUBEP3M


# sys.path.append(path_analy+'ToMesh')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
# # sys.path.insert(0, "/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh")
# import ToMeshCUBEP3M

# # filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# filepath = path_input
# # filepath = '/home/disraelcunha/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# redshift = '5.000'
# redshift_ = '5'
# Nmesh = 48
# BoxSize = 96
# nfiles = 8
# ncells = 96

# # mesh = ToMeshCUBEP3M.readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
# mesh = ToMeshCUBEP3M.Mesh_Wake(Nmesh,BoxSize)



# Input from Binary files

sys.path.append(path_analy+'read')
import Read_slices


# filepath = path_input
# redshift = '10'
# redshift_ = '10'
# Nmesh = 48
# BoxSize = 96
# nfiles = 8
# # nfiles = 512

# filepath = path_input
# redshift = '3'
# redshift_ = '3'
# Nmesh = 512
# BoxSize = 4
# nfiles = 32
# # nfiles = 512

filepath = path_input
redshift = '3'
redshift_ = '3'
Nmesh = 512
BoxSize = 4
nfiles = 1
# nfiles = 512

# mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)
mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,depth,filepath,redshift)



# # Plot Projection

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# sys.path.append(path_analy+'2d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# import Projection2d

# # save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+'.png'
# save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+anglid+'.png'
# # Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)
# proj_2d= Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)




# # Plot Projection slices

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
sys.path.append(path_analy+'2d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
import Projection2d

# slice_list = [0]
# # slice_list = list(range(32))
# sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]
# deept = 1

dept = 16
slice_list = list(range(int(depth/dept)))
sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]

save_plot_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]
test= Projection2d.plot_2d_proj_eachSlice(mesh,slice_list,dept,save_plot_2d_proj_fig_list)








# # histogram 3D

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append(path_analy+'3d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# import Analysis3d

# # save_3Dhist_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/3dhist_s1001_z0.png'
# # save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+'.png'
# save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+anglid+'.png'
# # save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+'.png'
# Analysis3d.histogram_3d(mesh,save_3Dhist_filename)





# # histogram 3D slices

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
sys.path.append(path_analy+'3d')
import Analysis3d

# # slice_list = [0]
# slice_list = list(range(32))
# sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]
# dept = 1

dept = 16
slice_list = list(range(int(depth/dept)))
sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]

save_3Dhist_2d_proj_fig_list = [ path_out+'3dhist_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]
Analysis3d.histogram_3d_eachSlice(mesh,slice_list,dept,save_3Dhist_2d_proj_fig_list)





# # Power Spectrum

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# sys.path.append(path_analy+'ps')
# import powerSpectrum_nbodykit

# # save_PS_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/PS_s1001_z0.png'
# save_PS_filename = path_out+'PS_'+sample+'_z'+redshift_+anglid+'.png'
# powerSpectrum_nbodykit.powerSpectrum(mesh,save_PS_filename)




# # compute the 2D power

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# sys.path.append(path_analy+'ps')
# import powerSpectrum_nbodykit

# nmu=5
# # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# save_PS2D_filename = path_out+'PS2D_'+sample+'_z'+redshift_+anglid+'.png'
# # Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
# Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu,save_PS2D_filename)



