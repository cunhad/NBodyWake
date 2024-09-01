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
        f = CustomDataFormatNbodykit.CUBEP3MCatalog(filename,Nmesh,ncells,nfiles)
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


def readCUBEP3M2(Nmesh,BoxSize,nfiles,ncells,filepath,redshift):
    
    import numpy as np
    from nbodykit.lab import ArrayMesh
    
    # generate random data on a 128^3 mesh
    
    Nmesh_x = Nmesh[0]
    Nmesh_y = Nmesh[1]
    Nmesh_z = Nmesh[2]
    
    grid = np.zeros((Nmesh_x,Nmesh_y,Nmesh_z))
    
    grid_spacing_x = ncells/Nmesh_x
    grid_spacing_y = ncells/Nmesh_y
    grid_spacing_z = ncells/Nmesh_z
    
    for i in range(0,nfiles):
        
        # print(i)
   
        filename = filepath+redshift+'xv'+str(i)+'.dat'
        node = int(i)
        nc = ncells
        number_node_dim = nfiles**(1./3)   
        k_node = np.floor(node/number_node_dim**2)
        res = np.floor(node % number_node_dim**2)
        j_node = np.floor(res/number_node_dim);
        i_node=res % number_node_dim
        
        data_xv = np.fromfile(filename, dtype=np.float32 , offset=4*(12)).reshape((-1,6))
        
        data_xv[:,0] = data_xv[:,0] + (nc/number_node_dim)*i_node
        data_xv[:,1] = data_xv[:,1] + (nc/number_node_dim)*j_node
        data_xv[:,2] = data_xv[:,2] + (nc/number_node_dim)*k_node
             
        for particle_position in data_xv[:, 0:3]:
            x, y, z = particle_position
        
            i1 = int(np.floor(x / grid_spacing_x)) 
            i2 = (i1 + 1) % Nmesh_x
            dx1 = x / grid_spacing_x - i1
            dx2 = 1 - dx1
            i1 = i1  % Nmesh_x

        
            j1 = int(np.floor(y / grid_spacing_y)) 
            j2 = (j1 + 1) % Nmesh_y
            dy1 = y / grid_spacing_y - j1
            dy2 = 1 - dy1
            j1 = j1  % Nmesh_y

        
            k1 = int(np.floor(z / grid_spacing_z)) 
            k2 = (k1 + 1) % Nmesh_z
            dz1 = z / grid_spacing_z - k1
            dz2 = 1 - dz1
            k1 = k1  % Nmesh_z

            
            # if (dx1 < 0 or dx2 < 0 or dy1 < 0 or dy2 < 0 or dz1 < 0 or dz2 < 0) and aux == 1:
            #     aux = 2
            #     print(x,y,z)
            #     print(dx1,dx2,dy1,dy2,dz1,dz2)
                
        
            # Update the grid with weights
            grid[i1, j1, k1] += dx1 * dy1 * dz1
            grid[i2, j1, k1] += dx2 * dy1 * dz1
            grid[i1, j2, k1] += dx1 * dy2 * dz1
            grid[i2, j2, k1] += dx2 * dy2 * dz1
            grid[i1, j1, k2] += dx1 * dy1 * dz2
            grid[i2, j1, k2] += dx2 * dy1 * dz2
            grid[i1, j2, k2] += dx1 * dy2 * dz2
            grid[i2, j2, k2] += dx2 * dy2 * dz2
        
    
    #normalize
    
    grid = grid/np.average(grid)
    
    # inititalize the mesh
    mesh = ArrayMesh(grid, BoxSize=BoxSize)
    
    return mesh.to_real_field(), grid, data_xv





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

    # mesh = RealField(mesh)
    
    return mesh.to_real_field()


# #Example Mesh_Wake

# from matplotlib import pyplot as plt

# Nmesh = 48
# BoxSize = 96
# mesh = Mesh_Wake(Nmesh,BoxSize)
# # plt.imshow(mesh.preview(axes=[0,1]))
# plt.imshow(mesh.preview(axes=[0,2]))






