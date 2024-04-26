#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 12:54:16 2023

@author: asus
"""


def read_slices_bin(Nmesh,BoxSize_,nfiles,filepath,redshift):
    
    from nbodykit.lab import ArrayMesh
    import numpy
    
    data = numpy.zeros(shape=(Nmesh,Nmesh,nfiles))
    
    for i in range(1,nfiles+1):
    # for i in range(1,1+1):
        filename = filepath+'data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/'+'_1_2dproj_z'+redshift+'_data_sl'+str(i)+'.bin'
        # print(filename)
        data_ = numpy.fromfile(filename, dtype=numpy.float32)
        data_ =   numpy.reshape(data_,(Nmesh,Nmesh)) 
        data[0:Nmesh,0:Nmesh,i-1] = data_

    # data =  data[0:32,0:32,0:32]  
    data = numpy.swapaxes(data, 1, 2)
    data = data/numpy.average(data)
    
    mesh = ArrayMesh(data, BoxSize=BoxSize_).to_real_field()
    
    
    return mesh
    # return data[0:Nmesh,0:Nmesh,0]
    # data_2 =   numpy.reshape(data_,(Nmesh,Nmesh)) 

    # return data_