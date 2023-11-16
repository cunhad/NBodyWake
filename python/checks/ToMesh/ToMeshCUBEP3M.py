#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 11:13:58 2023

@author: Disrael
"""


#from nbodykit.lab import cosmology

from nbodykit.lab import *
import CustomDataFormatNbodykit

# import numpy as np
from pmesh.pm import ParticleMesh, RealField







def readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift):
    
    
    # initialize the mesh
    
   
    pm = ParticleMesh(Nmesh=[Nmesh,Nmesh,Nmesh])
    mesh = RealField(pm)
    mesh[...] = 0.0
    
    for i in range(0,8):
    # for i in range(0,8):    
        filename = filepath+redshift+'xv'+str(i)+'.dat'
        # Reading a Custom Data Format¶
        f = CustomDataFormatNbodykit.CUBEP3MCatalog(filename,ncells,nfiles)
        #compute the \delta+1 values and substracts 1 to obtain the dc
        mesh = f.to_mesh(Nmesh,BoxSize).paint(mode='real')-1+mesh
        
    mesh = mesh/8
    
    return mesh+1


#Example readCUBEP3M


# from matplotlib import pyplot as plt

# filepath = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/'
# redshift = '0.000'
# Nmesh = 48
# BoxSize = 96
# nfiles = 8
# ncells = 96

# mesh = readCUBEP3M(Nmesh,BoxSize,nfiles,ncells,filepath,redshift)

# # plt.figure()    
# plt.imshow(mesh.preview(axes=[0,1]))








def Mesh_Wake(Nmesh,BoxSize_):
    
    from nbodykit.lab import ArrayMesh
    import numpy
    
    # generate random data on a 128^3 mesh
    data = numpy.random.random(size=(Nmesh,Nmesh,Nmesh))
    
    # print(int(numpy.floor(Nmesh/2)))
    
    for i in range(0,Nmesh):
        for j in range(0,Nmesh):
            data[i,j,int(Nmesh/2)]=data[i,j,int(Nmesh/2)]+1
    
    # inititalize the mesh
    mesh = ArrayMesh(data, BoxSize=BoxSize_)


    
    return mesh


# #Example Mesh_Wake

# from matplotlib import pyplot as plt

# Nmesh = 48
# BoxSize = 96
# mesh = Mesh_Wake(Nmesh,BoxSize)
# # plt.imshow(mesh.preview(axes=[0,1]))
# plt.imshow(mesh.preview(axes=[0,2]))






