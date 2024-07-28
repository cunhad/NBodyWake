#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:07:14 2023

@author: Disrael
"""


# # Input form CUBEP3M
# import sys


# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
# sys.path.insert(0, "/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh")
# import ToMeshCUBEP3M

# # filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# filepath = '/home/disraelcunha/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# redshift = '0.000'
# Nmesh = 48
# BoxSize = 96
# nfiles = 8
# ncells = 96

# # mesh = ToMeshCUBEP3M.readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
# mesh = ToMeshCUBEP3M.Mesh_Wake(Nmesh,BoxSize)



# Input from Binary files

import sys
import os



path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'read')
import Read_slices





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

path_input =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE8/4Mpc_2048c_1024p_zi63_nowakem/sample5046/data/1lf_0.5rf/NSIDE_8/anglid_1/-43-113--256pv_0.10211--0.62099-0.7854ra/2dproj/dm/"
path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE8/plots2/4Mpc_2048c_1024p_zi63_nowakem/"
sample = "sample5046"
anglid = "_aid1"

filepath = path_input
redshift = '3'
redshift_ = '3'
Nmesh = 512
BoxSize = 4
nfiles = 1
depth = 32

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

mesh= Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles, depth,filepath,redshift)




# # # Plot Projection

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# sys.path.append(path_analy+'2d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# import Projection2d

# save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+anglid+'.png'
# # Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)
# proj_2d= Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)




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






# # histogram 3D slices

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
sys.path.append(path_analy+'3d')
import Analysis3d

# # slice_list = [0]
# slice_list = list(range(32))
# sliceId =  ["_sl" + str(slice_list[i])  for i in slice_list]
# dept = 1

dept = 4
slice_list = list(range(int(32/dept)))
sliceId =  ["_sldp" + str(slice_list[i])  for i in slice_list]

save_3Dhist_2d_proj_fig_list = [ path_out+'2dproj_'+sample+'_z'+redshift_+anglid+sliceId[i]+'.png' for i in slice_list]
Analysis3d.histogram_3d_eachSlice(mesh,slice_list,dept,save_3Dhist_2d_proj_fig_list)








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
