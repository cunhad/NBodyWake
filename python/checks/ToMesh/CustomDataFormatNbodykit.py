#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 19:38:11 2023

@author: Disrael
"""

from nbodykit.io.base import FileType
import os
import numpy as np


class CUBEP3MFile(FileType):
    """
    A file-like object to read numpy ``.npy`` files
    """
    def __init__(self, path,nc,nnodes):
        self.path = path
        self.attrs = {}
        # load the data and set size and dtype
        self.nc = nc  #number of cells per dimension
        self.nnodes = nnodes  #number of nodes (each for file)
        # self._data = numpy.load(self.path)
        # self._data = np.fromfile(self.path,dtype=np.float32)
        self.size = int((os.path.getsize(path)/(4*8))-12) # total size
        # self.dtype = np.dtype(np.float32) # data dtype    
        self.dtype = np.dtype([('Position', ('f4', 3)),('Velocity', ('f4', 3))])
        # self.np = np #number of particles per dimension

    def read(self, columns, start, stop, step=1):
        """
        Read the specified column(s) over the given range
        """
        node = int(self.path.split('/')[-1].split('xv')[-1].split('.')[0])
        number_node_dim = self.nnodes**(1./3)
        
        k_node = np.floor(node/number_node_dim**2)
        res = np.floor(node % number_node_dim**2)
        j_node = np.floor(res/number_node_dim);
        i_node=res % number_node_dim
        
        
        
        data = np.fromfile(self.path, dtype=np.float32, count=(stop-start)*6 , offset=4*(12+start*6)).reshape((-1,6))

        if 'Position' in columns:
            aux=0
            data[:,0] = data[:,0] + (self.nc/number_node_dim)*i_node
            data[:,1] = data[:,1] + (self.nc/number_node_dim)*j_node
            data[:,2] = data[:,2] + (self.nc/number_node_dim)*k_node
            


        if 'Velocity' in columns:
            aux=1
        
        
        # return data[start:stop:step,0:3]
        return data[:,0+3*aux:3+4*aux]
    
# Reading a Custom Data Format¶

from nbodykit.source.catalog.file import FileCatalogFactory

CUBEP3MCatalog = FileCatalogFactory('CUBEP3MCatalog', CUBEP3MFile)


# # Reading a Custom Data Format¶


# f = CUBEP3MCatalog('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv7.dat',96,8)

    

# print(f)
# print("columns = ", f.columns) # default Weight,Selection also present
# print("total size = ", f.csize)


# print("pos = ", f['Position'])
# pos_b =  f['Position']
# pos_bb = pos_b.compute()

# print("vel_b = ", f['Velocity'])
# vel_b =  f['Velocity']
# vel_bb = vel_b.compute()

# print(f.size)
