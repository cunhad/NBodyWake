#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 12:54:16 2023

@author: asus
"""


def read_slices_bin(Nmesh,BoxSize_,nfiles,filepath,redshift):
    
    from nbodykit.lab import ArrayMesh
    import numpy
    
    Nmesh_x = Nmesh[0]
    Nmesh_y = Nmesh[1]
    Nmesh_z = Nmesh[2]
    
    
    # data = numpy.zeros(shape=(Nmesh,Nmesh,depth))
    
    if nfiles==1:        
        filename = filepath+'_1_2dproj_z'+redshift+'_data_slAll.bin'
        # filename = filepath+'_6_2dproj_z'+redshift+'_data_slAll.bin'
        data = numpy.fromfile(filename, dtype=numpy.float32)
        data =   numpy.reshape(data,(Nmesh_x,Nmesh_y,Nmesh_z)) 
        # data = numpy.swapaxes(data, 0,1)
        data = numpy.swapaxes(data, 1,2)
    else:
        for i in range(1,nfiles+1):
        # for i in range(1,1+1):
            filename = filepath+'_1_2dproj_z'+redshift+'_data_sl'+str(i)+'.bin'
            # print(filename)
            data_ = numpy.fromfile(filename, dtype=numpy.float32)
            data_ =   numpy.reshape(data_,(Nmesh_y,Nmesh_z)) 
            data[i-1,0:Nmesh_z,0:Nmesh_z] = data_
        data = numpy.swapaxes(data, 1, 2)
    
    data = data/numpy.average(data)    
    mesh = ArrayMesh(data, BoxSize=BoxSize_).to_real_field()
    
    
    return mesh
    # return data[0:Nmesh,0:Nmesh,0]
    # data_2 =   numpy.reshape(data_,(Nmesh,Nmesh)) 

    # return data_
    # return data    
    

def extract_numbers(input_string):
    
    import re

    
    # Find all patterns that start with "-" or "--" followed by a number, or just a number
    matches = re.findall(r'(?:(?:-{1,2})?\d+\.?\d*)', input_string)

    # Convert the matches to the correct sign based on the length of the "-" symbols
    results = []
    for match in matches:
        if match.startswith('--'):
            results.append(-float(match[2:]))  # Extract number after "--" and make it negative
        elif match.startswith('-'):
            results.append(float(match[1:]))  # Extract number after "-" and keep it positive
        else:
            results.append(float(match))      # No "-" means it's positive
    
    return results

def extract_pv_ra(string):
    
    import re

    
    # Extract the part before 'pv' and the part between 'pv' and 'ra'
    
    
    
    string = string.split('/')[-4]

    
    string = '_'.join(string.split("_")[-2:])
    
    pv_part = re.search(r'(.*)pv', string).group(1)
    ra_part = re.search(r'pv_(.*)ra', string).group(1)

    # Extract the numbers for pv, handling '--' as a negative sign only for the appropriate part
    pv = extract_numbers(pv_part) 

    # Extract the numbers for ra
    ra = extract_numbers(ra_part) 


    return pv, ra    