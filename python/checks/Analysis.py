#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:07:14 2023

@author: Disrael
"""


# import ToMesh/ToMesh

import sys


# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh')
# sys.path.insert(0, "/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ToMesh")
import ToMeshCUBEP3M

# filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
filepath = '/home/disraelcunha/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
redshift = '0.000'
Nmesh = 48
BoxSize = 96
nfiles = 8
ncells = 96

# mesh = ToMeshCUBEP3M.readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)
mesh = ToMeshCUBEP3M.Mesh_Wake(Nmesh,BoxSize)


# # Plot Projection

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/2d')
# import Projection2d

# save_plot_2d_proj_fig = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/2dproj_s1001_z0.png'
# Projection2d.plot_2d_proj(mesh,save_plot_2d_proj_fig)




# # histogram 3D

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# # sys.path.append('/home/disraelcunha/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/3d')
# import Analysis3d

# save_3Dhist_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/3dhist_s1001_z0.png'
# Analysis3d.histogram_3d(mesh,save_3Dhist_filename)




# # Power Spectrum

# sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
# import powerSpectrum_nbodykit

# save_PS_filename = '/home/asus/Dropbox/extras/storage/graham/small_res/plots/64Mpc_96c_48p_zi255_nowakem/PS_s1001_z0.png'
# powerSpectrum_nbodykit.powerSpectrum(mesh,save_PS_filename)




# compute the 2D power

sys.path.append('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
import powerSpectrum_nbodykit

nmu=5
# Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)
Pkmu = powerSpectrum_nbodykit.powerSpectrum2d(mesh,nmu)



