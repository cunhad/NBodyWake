#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 12:54:16 2023

@author: asus
"""


def read_slices_bin(Nmesh,BoxSize_,nfiles,filepath,redshift):
    
    # from nbodykit.lab import ArrayMesh
    import numpy
    import glob

    
    Nmesh_x = int(Nmesh[0])
    Nmesh_y = int(Nmesh[1])
    Nmesh_z = int(Nmesh[2])
    
    
    # data = numpy.zeros(shape=(Nmesh,Nmesh,depth))
    
    if nfiles==1:    
        # Define the file pattern with a wildcard
        file_pattern = filepath + '_*_2dproj_z' + redshift + '_data_slAll.bin'
        file_list = glob.glob(file_pattern)
        filename = file_list[0] 
        # print(filename)
        
        # filename = filepath+'_1_2dproj_z'+redshift+'_data_slAll.bin'
        # filename = filepath+'_6_2dproj_z'+redshift+'_data_slAll.bin'
        # filename = filepath+'__2dproj_z'+redshift+'_data_slAll.bin'
        
        data = numpy.fromfile(filename, dtype=numpy.float32)
        data =   numpy.reshape(data,(Nmesh_z,Nmesh_x,Nmesh_y,)) 
        # data = numpy.swapaxes(data, 1,2)
        # data = numpy.swapaxes(data, 0,1)
        # data = numpy.swapaxes(data, 1,2)
        # # data = numpy.transpose(data, (1, 2, 0))
        data = numpy.moveaxis(data, [0, 1, 2], [2, 1, 0])
    else:
        # have to fix this!
        data = numpy.zeros((nfiles, Nmesh_z, Nmesh_z))
        for i in range(1,nfiles+1):
        # for i in range(1,1+1):
            filename = filepath+'_1_2dproj_z'+redshift+'_data_sl'+str(i)+'.bin'
            # print(filename)
            data_ = numpy.fromfile(filename, dtype=numpy.float32)
            data_ =   numpy.reshape(data_,(Nmesh_y,Nmesh_z)) 
            data[i-1,0:Nmesh_y,0:Nmesh_z] = data_
        data = numpy.swapaxes(data, 1, 2)
    
    data = (data/numpy.average(data)) - 1
    # mesh = ArrayMesh(data, BoxSize=BoxSize_).to_real_field()
    
    
    # return mesh
    # return data[0:Nmesh,0:Nmesh,0]
    # data_2 =   numpy.reshape(data_,(Nmesh,Nmesh)) 

    # return data_
    return data    
    

# def extract_numbers(input_string):
    
#     import re

    
#     # Find all patterns that start with "-" or "--" followed by a number, or just a number
#     matches = re.findall(r'(?:(?:-{1,2})?\d+\.?\d*)', input_string)

#     # Convert the matches to the correct sign based on the length of the "-" symbols
#     results = []
#     for match in matches:
#         if match.startswith('--'):
#             results.append(-float(match[2:]))  # Extract number after "--" and make it negative
#         elif match.startswith('-'):
#             results.append(float(match[1:]))  # Extract number after "-" and keep it positive
#         else:
#             results.append(float(match))      # No "-" means it's positive
    
#     return results

# def extract_pv_ra(string):
    
#     import re

    
#     # Extract the part before 'pv' and the part between 'pv' and 'ra'
    
    
    
#     string = string.split('/')[-4]

    
#     string = '_'.join(string.split("_")[-2:])
    
#     pv_part = re.search(r'(.*)pv', string).group(1)
#     ra_part = re.search(r'pv_(.*)ra', string).group(1)

#     # Extract the numbers for pv, handling '--' as a negative sign only for the appropriate part
#     pv = extract_numbers(pv_part) 

#     # Extract the numbers for ra
#     ra = extract_numbers(ra_part) 


#     return pv, ra    

def extract_numbers(input_string, is_first_number=False):
    
    import re
    
    # Find all patterns that start with "-" or "--" followed by a number, or just a number
    matches = re.findall(r'(?:(?:-{1,2})?\d+\.?\d*)', input_string)

    # Convert the matches to the correct sign based on the length of the "-" symbols
    results = []
    for i, match in enumerate(matches):
        if i == 0 and is_first_number:
            # First number logic: "-" means negative, nothing means positive
            if match.startswith('-'):
                results.append(-float(match[1:]))  # Negative number
            else:
                results.append(float(match))       # Positive number
        else:
            # Other numbers logic: "--" means negative, "-" means positive
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
    string = string.split('/')[-4]  # Adjusted to correctly get the pattern from filepath
    
    string = '_'.join(string.split("_")[-2:])
    
    pv_part = re.search(r'(.*)pv', string).group(1)
    ra_part = re.search(r'pv_(.*)ra', string).group(1)
    
    # pv_part = re.search(r'(.*)pv', string).group(1)
    # ra_part = re.search(r'pv_(.*)ra', string).group(1)

    # Extract the numbers for pv, with special logic for the first number
    pv = extract_numbers(pv_part, is_first_number=True) 

    # Extract the numbers for ra, with special logic for the first number
    ra = extract_numbers(ra_part, is_first_number=True)

    return pv, ra


def filepath_from_wake_to_nowake(filepath):

    import re


    # Replace the pattern "*Mpc_*c_*p_zi*_.*" with "*Mpc_*c_*p_zi*_nowakem"
    filepath = re.sub(r'(\d+Mpc_\d+c_\d+p_zi\d+)_\w+', r'\1_nowakem', filepath)
    
    # Remove the occurrence of 'half_lin_cutoff_half_tot_pert_nvpw*'
    # filepath = filepath.replace('half_lin_cutoff_half_tot_pert_nvpw/', '')
    filepath = re.sub(r'half_lin_cutoff_half_tot_pert_nvpw\w*', '', filepath)

    
    return filepath

# def read_slice_bin(Nmesh,BoxSize_,filepath,redshift,slice_idx):
    
#     # from nbodykit.lab import ArrayMesh
#     import numpy as np
#     import glob

    
#     Nmesh_z, Nmesh_y, Nmesh_x = Nmesh
    
#     # elements_per_slice = Nmesh_y * Nmesh_z
#     # Calculate the offset in bytes for the desired slice
#     # Since each float32 is 4 bytes, multiply by 4 to get the byte offset
#     # offset = slice_idx * elements_per_slice * 4  # offset in bytes
   
#     # Define the file pattern with a wildcard
#     file_pattern = filepath + '_*_2dproj_z' + redshift + '_data_slAll.bin'
#     file_list = glob.glob(file_pattern)
#     filename = file_list[0] 
    
#     # filename = filepath+'_1_2dproj_z'+redshift+'_data_slAll.bin'
#     # filename = filepath+'_6_2dproj_z'+redshift+'_data_slAll.bin'
#     # filename = filepath+'__2dproj_z'+redshift+'_data_slAll.bin'
    
#     # Initialize an empty array to store the 2D slice
#     slice_2d = np.zeros((Nmesh_y, Nmesh_z), dtype=np.float32)
    
#     # Each row has Nmesh_x elements, so each slice is separated by Nmesh_y * Nmesh_x elements
#     # row_size_bytes = Nmesh_x * 4  # each element is 4 bytes (float32)
#     # slice_stride = Nmesh_x * Nmesh_y * Nmesh_z * 4  # full 3D array size in bytes for a single slice
    
#     # Open the file and read each row of the desired slice
#     with open(filename, 'rb') as file:
#         for j in range(Nmesh_y):
#             for k in range(Nmesh_z):
#                 offset = (slice_idx + ((k * Nmesh_y) + j) * Nmesh_x)*4
#                 # Calculate the starting byte offset for this row in the specific slice
            
#                 file.seek(offset)
            
#                 # Read the row and store it in the corresponding location in the 2D slice
#                 slice_2d[j, k] = np.fromfile(file, dtype=np.float32, count=1)
    
    
#     # data = numpy.fromfile(filename, dtype=numpy.float32)
#     # data =   numpy.reshape(data,(Nmesh_z,Nmesh_x,Nmesh_y,)) 
#     # # data = numpy.swapaxes(data, 1,2)
#     # # data = numpy.swapaxes(data, 0,1)
#     # # data = numpy.swapaxes(data, 1,2)
#     # # # data = numpy.transpose(data, (1, 2, 0))
#     # data = numpy.moveaxis(data, [0, 1, 2], [2, 1, 0])
   




    
#     # data = (data/numpy.average(data)) - 1
#     # # mesh = ArrayMesh(data, BoxSize=BoxSize_).to_real_field()
    
    
#     # return mesh
#     # return data[0:Nmesh,0:Nmesh,0]
#     # data_2 =   numpy.reshape(data_,(Nmesh,Nmesh)) 

#     # return data_
#     return slice_2d   


def read_slice_bin(Nmesh,BoxSize_,filepath,redshift,slice_idx):
    
    # from nbodykit.lab import ArrayMesh
    import numpy as np
    import glob

    Nmesh_x, Nmesh_y, Nmesh_z = Nmesh
    
    # Calculate the number of elements in a 2D slice (Nmesh_x * Nmesh_y)
    elements_per_slice = Nmesh_x * Nmesh_y
    
    # Calculate the offset in bytes for the desired slice
    # Since each float32 is 4 bytes, multiply by 4 to get the byte offset
    offset = slice_idx * elements_per_slice * 4  # offset in bytes
    
    
    
   
   
    # Define the file pattern with a wildcard
    file_pattern = filepath + '_*_2dproj_z' + redshift + '_data_slAll.bin'
    file_list = glob.glob(file_pattern)
    filename = file_list[0] 
    
    # Open the file and read only the necessary slice
    with open(filename, 'rb') as file:
        # Move the file pointer to the calculated offset
        file.seek(offset)
        
        # Read only the required number of elements for one slice
        slice_2d = np.fromfile(file, dtype=np.float32, count=elements_per_slice)
        # Reshape to 2D array (Nmesh_x, Nmesh_y)
        
    slice_2d = slice_2d.reshape((Nmesh_x, Nmesh_y))
    slice_2d = slice_2d.T
    

    return slice_2d   