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
    grid = np.zeros((Nmesh,Nmesh,Nmesh))
    
    grid_spacing = ncells/Nmesh
    
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
        
        # for particle_position in data_xv[:,0:3]:
        #     x, y, z = particle_position
        #     # cell_index_x = int(x / grid_spacing) % Nmesh
        #     # cell_index_y = int(y / grid_spacing) % Nmesh
        #     # cell_index_z = int(z / grid_spacing) % Nmesh
            
        #     # cell_index_x = int(np.floor(x / grid_spacing)) % Nmesh
        #     # cell_index_y = int(np.floor(y / grid_spacing)) % Nmesh
        #     # cell_index_z = int(np.floor(z / grid_spacing)) % Nmesh

        #     cell_index_x = int(np.floor(x / grid_spacing)) 
        #     cell_index_y = int(np.floor(y / grid_spacing)) 
        #     cell_index_z = int(np.floor(z / grid_spacing)) 
            
        #     # Calculate fractional contributions to neighboring cells
        #     dx = x / grid_spacing - cell_index_x
        #     dy = y / grid_spacing - cell_index_y
        #     dz = z / grid_spacing - cell_index_z    
            
        #     for i in range(cell_index_x, cell_index_x + 2):
        #         i = i % Nmesh
        #         for j in range(cell_index_y, cell_index_y + 2):
        #             j = j % Nmesh
        #             for k in range(cell_index_z, cell_index_z + 2):
        #                 k = k % Nmesh
                    
        #                 # weight = (1 - dx) * (1 - dy) * (1 - dz) if i == cell_index_x else dx * (1 - dy) * (1 - dz)
        #                 # weight += (1 - dx) * dy * (1 - dz) if j == cell_index_y else dx * dy * (1 - dz)
        #                 # weight += (1 - dx) * (1 - dy) * dz if k == cell_index_z else dx * dy * dz
        #                 # grid[i, j, k] += weight
        #                 if 0 <= i < Nmesh and 0 <= j < Nmesh and 0 <= k < Nmesh:
        #                       weight = (1 - dx) * (1 - dy) * (1 - dz) if i == cell_index_x else dx * (1 - dy) * (1 - dz)
        #                       weight += (1 - dx) * dy * (1 - dz) if j == cell_index_y else dx * dy * (1 - dz)
        #                       weight += (1 - dx) * (1 - dy) * dz if k == cell_index_z else dx * dy * dz
        #                       grid[i, j, k] += weight
        
        # for particle_position in data_xv[:, 0:3]:
        #     x, y, z = particle_position
        
        #     # Find the indices of the cell lower corner
        #     cell_index_x = int(np.floor(x / grid_spacing))
        #     cell_index_y = int(np.floor(y / grid_spacing))
        #     cell_index_z = int(np.floor(z / grid_spacing))
            
        #     # Calculate fractional contributions to neighboring cells
        #     dx = x / grid_spacing - cell_index_x
        #     dy = y / grid_spacing - cell_index_y
        #     dz = z / grid_spacing - cell_index_z
            
        #     # Iterate over the 8 corners of the cell
        #     for i in range(2):
        #         for j in range(2):
        #             for k in range(2):
        #                 # Adjust indices for the 8 neighboring cells
        #                 ix = (cell_index_x + i) % Nmesh
        #                 iy = (cell_index_y + j) % Nmesh
        #                 iz = (cell_index_z + k) % Nmesh
        
        #                 # Calculate weight for each corner based on the relative position
        #                 weight = ((1 - i) * (1 - dx) + i * dx) * \
        #                          ((1 - j) * (1 - dy) + j * dy) * \
        #                          ((1 - k) * (1 - dz) + k * dz)
                        
        #                 # Add weight to the appropriate grid cell
        #                 if 0 <= ix < Nmesh and 0 <= iy < Nmesh and 0 <= iz < Nmesh:
        #                     grid[ix, iy, iz] += weight        
                
        # for particle_position in data_xv[:,0:3]:
        #     x, y, z = particle_position
        #     # i1 = int(np.floor(x / grid_spacing - 0.5) % Nmesh)
        #     # i2 = i1 + 1
        #     # dx1=(i1-x)
        #     # dx2=1-dx1
            
        #     i1 = int(np.floor(x / grid_spacing - 0.5) % Nmesh)
        #     i2 = i1 + 1
        #     dx1=(i1-x)
        #     dx2=1-dx1
            
        #     j1 = int(np.floor(y / grid_spacing - 0.5) % Nmesh)
        #     j2 = j1 + 1
        #     dy1=j1-y
        #     dy2=1-dy1        
            
        #     k1 = int(np.floor(z / grid_spacing - 0.5) % Nmesh)
        #     k2 = k1 + 1
        #     dz1=k1-z
        #     dz2=1-dz1
            
        #     if 0 <= i1 < Nmesh and 0 <= j1 < Nmesh and 0 <= k1 < Nmesh and 0 <= i2 < Nmesh and 0 <= j2 < Nmesh and 0 <= k2 < Nmesh:
                
        #         grid[i1,j1,k1]=grid[i1,j1,k1]+dx1*dy1*dz1
        #         grid[i2,j1,k1]=grid[i2,j1,k1]+dx2*dy1*dz1
        #         grid[i1,j2,k1]=grid[i1,j2,k1]+dx1*dy2*dz1
        #         grid[i2,j2,k1]=grid[i2,j2,k1]+dx2*dy2*dz1
        #         grid[i1,j1,k2]=grid[i1,j1,k2]+dx1*dy1*dz2
        #         grid[i2,j1,k2]=grid[i2,j1,k2]+dx2*dy1*dz2
        #         grid[i1,j2,k2]=grid[i1,j2,k2]+dx1*dy2*dz2
        #         grid[i2,j2,k2]=grid[i2,j2,k2]+dx2*dy2*dz2
        
        # for particle_position in data_xv[:,0:3]:
        #     x, y, z = particle_position
        #     # i1 = int(np.floor(x / grid_spacing - 0.5) % Nmesh)
        #     # i2 = i1 + 1
        #     # dx1=(i1-x)
        #     # dx2=1-dx1
            
        #     i1 = int(np.floor(x / grid_spacing - 0.5) + 1 ) % Nmesh
        #     i2 = (i1 + 1) % Nmesh
        #     dx1 = i1- (x / grid_spacing)
        #     dx2 = 1-dx1
            
        #     j1 = int(np.floor(y / grid_spacing - 0.5) +1 ) % Nmesh
        #     j2 = (j1 + 1) % Nmesh
        #     dy1 = j1-(y / grid_spacing)
        #     dy2 = 1-dy1        
            
        #     k1 = int(np.floor(z / grid_spacing - 0.5) +1)  % Nmesh
        #     k2 = (k1 + 1) % Nmesh
        #     dz1 = k1-(z / grid_spacing)
        #     dz2 = 1-dz1
            
        #     # if 0 <= i1 < Nmesh and 0 <= j1 < Nmesh and 0 <= k1 < Nmesh and 0 <= i2 < Nmesh and 0 <= j2 < Nmesh and 0 <= k2 < Nmesh:
                
        #     grid[i1,j1,k1]=grid[i1,j1,k1]+dx1*dy1*dz1
        #     grid[i2,j1,k1]=grid[i2,j1,k1]+dx2*dy1*dz1
        #     grid[i1,j2,k1]=grid[i1,j2,k1]+dx1*dy2*dz1
        #     grid[i2,j2,k1]=grid[i2,j2,k1]+dx2*dy2*dz1
        #     grid[i1,j1,k2]=grid[i1,j1,k2]+dx1*dy1*dz2
        #     grid[i2,j1,k2]=grid[i2,j1,k2]+dx2*dy1*dz2
        #     grid[i1,j2,k2]=grid[i1,j2,k2]+dx1*dy2*dz2
        #     grid[i2,j2,k2]=grid[i2,j2,k2]+dx2*dy2*dz2
        # aux = 1
        for particle_position in data_xv[:, 0:3]:
            x, y, z = particle_position
        
            i1 = int(np.floor(x / grid_spacing)) 
            i2 = (i1 + 1) % Nmesh
            dx1 = x / grid_spacing - i1
            dx2 = 1 - dx1
            i1 = i1  % Nmesh

        
            j1 = int(np.floor(y / grid_spacing)) 
            j2 = (j1 + 1) % Nmesh
            dy1 = y / grid_spacing - j1
            dy2 = 1 - dy1
            j1 = j1  % Nmesh

        
            k1 = int(np.floor(z / grid_spacing)) 
            k2 = (k1 + 1) % Nmesh
            dz1 = z / grid_spacing - k1
            dz2 = 1 - dz1
            k1 = k1  % Nmesh

            
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






