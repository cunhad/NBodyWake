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
    # /home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/ /home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/
    


# splited = sys.argv[0].split('/')   
# path_analy = "/".join(splited[0:-1])+'/'
# # path_analy =  sys.argv[1]
path_analy = os.getcwd() +'/' 
path_input =  sys.argv[1]
path_out = sys.argv[2]

if sys.argv[1].split('/')[-2][0] == 's':
    sample = sys.argv[1].split('/')[-2]
else:
    sample = sys.argv[1].split('/')[-3]
    


sys.path.append(path_analy+'ToMesh')
# sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
# sys.path.insert(0, "/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh")
import ToMeshCUBEP3M

# filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
filepath = path_input
# filepath = '/home/disraelcunha/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
redshift = '5.000'
redshift_ = '5'
Nmesh = 48
BoxSize = 96
nfiles = 8
ncells = 96

mesh = ToMeshCUBEP3M.readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
# mesh = ToMeshCUBEP3M.Mesh_Wake(Nmesh,BoxSize)


# # Plot Projection

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# sys.path.append(path_analy+'2d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# import Projection2d

# save_plot_2d_proj_fig = path_out+'2dproj_'+sample+'_z'+redshift_+'.png'
# # Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)
# proj_2d= Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)




# # histogram 3D

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# sys.path.append(path_analy+'3d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# import Analysis3d

# # save_3Dhist_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/3dhist_s1001_z0.png'
# save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+'.png'

# # save_3Dhist_filename = path_out+'3dhist_'+sample+'_z'+redshift_+'.png'
# Analysis3d.histogram_3d(mesh,save_3Dhist_filename)




# # Power Spectrum

# # sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# sys.path.append(path_analy+'ps')
# import powerSpectrum_nbodykit

# # save_PS_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/PS_s1001_z0.png'
# save_PS_filename = path_out+'PS_'+sample+'_z'+redshift_+'.png'
# powerSpectrum_nbodykit.powerSpectrum(mesh,save_PS_filename)




# compute the 2D power

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
sys.path.append(path_analy+'ps')
import powerSpectrum_nbodykit

nmu=5
# Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
save_PS2D_filename = path_out+'PS2D_'+sample+'_z'+redshift_+'.png'
# Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu,save_PS2D_filename)



