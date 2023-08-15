#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:15:03 2023

@author: asus
"""




#from nbodykit.lab import cosmology

from nbodykit.lab import *




import numpy as np
from matplotlib import pyplot as plt


import CustomDataFormatNbodykit


# Reading a Custom Data Format¶

Nmesh = 48
BoxSize = 96
nfiles = 8
ncells = 96

f = CustomDataFormatNbodykit.CUBEP3MCatalog('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat',ncells,nfiles)

    

# mesh = f.to_mesh(Nmesh=64,BoxSize=64).paint(mode='real')
# # plt.figure()    
# plt.imshow(mesh.preview(axes=[0,1]))

# # rfield_ = mesh.preview()

# mesh2 = mesh+1

# # plt.figure()    
# plt.imshow(mesh2.preview(axes=[0,1]))

# a1 = mesh.preview().max()
# a2 = mesh2.preview().max()
# pos_b =  src['Position']
# pos_bb = pos_b.compute()



# initialize the mesh

from pmesh.pm import ParticleMesh, RealField, ComplexField

pm = ParticleMesh(Nmesh=[Nmesh,Nmesh,Nmesh])
mesh = RealField(pm)
mesh[...] = 0.0

# mesh_ = mesh.preview()
# mesh2 = mesh_

# plt.imshow(mesh.preview(axes=[0,1]))


for i in range(0,8):
# for i in range(0,8):    
    filename = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv'+str(i)+'.dat'
    f = CustomDataFormatNbodykit.CUBEP3MCatalog(filename,ncells,nfiles)
    #compute the \delta+1 values and substracts 1 to obtain the dc
    mesh = f.to_mesh(Nmesh,BoxSize).paint(mode='real')-1+mesh
    
mesh = mesh/8
    
# plt.figure()    
plt.imshow(mesh.preview(axes=[0,1]))

a = mesh.preview().min()






