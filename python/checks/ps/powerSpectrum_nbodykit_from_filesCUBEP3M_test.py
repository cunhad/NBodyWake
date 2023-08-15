#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:15:03 2023

@author: asus
"""




#from nbodykit.lab import cosmology

from nbodykit.lab import *




#first read the files from GADGET


# # Plaintext Data

# import numpy
# from nbodykit.source.catalog import CSVCatalog

# # generate some fake ASCII data
# data = numpy.random.random(size=(100,5))

# # save to a plaintext file
# numpy.savetxt('/home/asus/Dropbox/extras/storage/test/catalog/csv-example.txt', data, fmt='%.7e')

# # name each of the 5 input columns
# names =['x', 'y', 'z', 'w', 'v']

# # read the data
# f = CSVCatalog('/home/asus/Dropbox/extras/storage/test/catalog/csv-example.txt', names)



# # combine x, y, z to Position, and add boxsize
# f['Position'] = f['x'][:, None] * [1, 0, 0] + f['y'][:, None] * [0, 1, 0] + f['z'][:, None] * [0, 0, 1]
# f.attrs['BoxSize'] = 1.0

# print(f)
# print("columns = ", f.columns) # default Weight,Selection also present
# print("total size = ", f.csize)



# #  read Binary Data (format understandable to NBodyKit)

# from nbodykit.source.catalog import BinaryCatalog

# # generate some fake data and save to a binary file
# with open('/home/asus/Dropbox/extras/storage/test/catalog/binary-example.dat', 'wb') as ff:
#     pos = numpy.random.random(size=(30, 3)) # fake Position column
#     vel = numpy.random.random(size=(30, 3)) # fake Velocity column
#     pos.tofile(ff); vel.tofile(ff); ff.seek(0)

# # create the binary catalog
# f = BinaryCatalog('/home/asus/Dropbox/extras/storage/test/catalog/binary-example.dat', [('Position', ('f8', 3)), ('Velocity', ('f8', 3))], size=30)

# print(f)
# print("columns = ", f.columns) # default Weight,Selection also present
# print("total size = ", f.csize)

# print("pos = ", f['Position'])
# pos_b =  f['Position']






# # read built in binary (from read example)

# import numpy as np


# data = np.fromfile('/home/asus/Dropbox/extras/storage/test/catalog/binary-example.dat', dtype=np.float64)


# # print(np.float32)







# # read built in binary  (from read CUBEP3M)

# import numpy as np

# data = np.fromfile('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat', dtype=np.float32)

# data = data[12:]        #skip header

# numParticles = int(len(data)/6)     #each particle has 6 numbers (x,y,z,vx,vy,vz)

# pos = np.ndarray(shape=(numParticles, 3), dtype=np.float32) 
# vel = np.ndarray(shape=(numParticles, 3), dtype=np.float32) 

# pos[:,0]=data[0::6]
# pos[:,1]=data[1::6]
# pos[:,2]=data[2::6]


# vel[:,0]=data[3::6]
# vel[:,1]=data[4::6]
# vel[:,2]=data[5::6]





# #  read Binary Data (from read CUBEP3M)



# create the binary catalog
# f = BinaryCatalog('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat', [('Position', ('f4', 3)), ('Velocity', ('f4', 3))], header_size=12*4)


# import os
# size = int((os.path.getsize('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat')/32)-12)

# import numpy as np
# data=np.dtype(np.float32)
# # data=np.float32
# data2 = data.names


# print(f)
# print("columns = ", f.columns) # default Weight,Selection also present
# print("total size = ", f.csize)

# print("pos_b = ", f['Position'])
# pos_b =  f['Position']
# pos_bb = pos_b.compute()

# print("vel_b = ", f['Velocity'])
# vel_b =  f['Velocity']
# vel_bb = vel_b.compute()


# import numpy as np

# print(np.float32)



# test 

import numpy as np
# # data[:,1]=np.linspace(1,15395,15395)



# data = np.ndarray(shape=(15395, 6), dtype=np.float32)
# data[:,0]=0
# data[:,1]=1
# data[:,2]=2
# data[:,3]=3
# data[:,4]=4
# data[:,5]=5


# data2=data[:,0:3]
# start = 0;
# stop=15395
# aux=1

# count_ = (stop-start)*6        
# data_in = np.fromfile('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat', dtype=np.float32, count=count_, offset=4*(12+start*6))
# # data1 = data_in[0::6]
# # data2 = data_in[1::6]
# # data3 = data_in[2::6]
# data = data_in.reshape((-1,6))



# data = np.fromfile('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat', dtype=np.float32, count=count_, offset=4*(12+start*6)).reshape((-1,6))
# data_in = data[:,0+3*aux:3+4*aux]

# f = CUBEP3MCatalog('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat')



# filename = '/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat'
# rank = int(filename.split('/')[-1].split('xv')[-1].split('.')[0])

# node = 6
# number_node_dim = 2
# k_node = np.floor(node/number_node_dim**2)
# res = np.floor(node % number_node_dim**2)
# j_node = np.floor(res/number_node_dim);
# i_node=res % number_node_dim



# to use CUBEP3MCatalog

# import sys
# sys.path.insert(1, '/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps')
import CustomDataFormatNbodykit

from matplotlib import pyplot as plt


# from /home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/python/checks/ps/ import CustomDataFormatNbodykit.py

# from nbodykit.source.catalog.file import FileCatalogFactory

# CUBEP3MCatalog = CustomDataFormatNbodykit.FileCatalogFactory('CUBEP3MCatalog', CUBEP3MFile)


# Reading a Custom Data Format¶

Nmesh = 48
BoxSize = 96
nfiles = 8
ncells = 96

f = CustomDataFormatNbodykit.CUBEP3MCatalog('/home/asus/Dropbox/extras/storage/graham/small_res/64Mpc_96c_48p_zi255_nowakem/sample1001/0.000xv0.dat',ncells,nfiles)



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

mesh = f.to_mesh(Nmesh,BoxSize).paint(mode='real')

# plt.figure()    
plt.imshow(mesh.preview(axes=[0,1]))


# data = numpy.ones(1500, dtype=[
#   ('Position', ('f4', 3)),
#   ('Velocity', ('f4', 3))]
#   )
# # data['Position'] = 0*data['Position'] 
# src = ArrayCatalog(data)
# mesh = src.to_mesh(Nmesh=2,BoxSize=2)       #computes the \delta+1 values
# rfield_ = mesh.preview()
# rfield = mesh.paint(mode='real').preview()


# a = np.sum(np.sum(np.sum(rfield_,axis = 0),axis = 0),axis = 0)  









