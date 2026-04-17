#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:31:41 2024

@author: asus
"""



# parse stuff

import argparse
import ast



def parse_range(value):
    if not value or not value.strip():
        return []
    
    numbers = set()
    parts = value.split(',')

    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            numbers.update(range(start, end + 1))
        else:
            numbers.add(int(part))

    return sorted(numbers)


parser = argparse.ArgumentParser(description='void anlysis')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--folderpath', type=str, help='')
parser.add_argument('--sample', type=str, help='')
parser.add_argument('--sampleC', type=str, default='', help='')
parser.add_argument('--path_out', type=str, help='')
parser.add_argument('--redshift', type=str, default='3', help='')
parser.add_argument('--Nmesh', default='[512,512,32]', help='')
parser.add_argument('--BoxSize', type=float, default=4, help='')
parser.add_argument('--nfiles', type=int, default=1, help='')
parser.add_argument('--resol_factor', type=float, default=0.5, help='')
parser.add_argument('--Nside', type=int,default=8, help='')
parser.add_argument('--rangeAng', type=parse_range, default='1-96', help="Range of angular ids (e.g., '15-20,50,80')")
parser.add_argument('--downSample', type=int, default=1, help="the donwsample cell number to sum up (e.g. ==2**k, k int)")



# parser.add_argument('--num_epochs', type=int, default=10, help='')

args = parser.parse_args()



# parameters

folderpath = args.folderpath
print("File Path in = "+ str(folderpath))

sample = args.sample
print("sample = "+ str(sample))

#aux sample name (eg.wake prensence characterization)
sampleC = args.sampleC
print("sampleC = "+ str(sampleC))

path_out = args.path_out
print("Path out = "+ str(path_out))

redshift = args.redshift
print("redshift= "+ str(redshift))

Nmesh =  ast.literal_eval(args.Nmesh)
print("Nmesh= "+ str(Nmesh))

BoxSize = args.BoxSize
print("BoxSize= "+ str(BoxSize))

nfiles = args.nfiles
print("nfiles= "+ str(nfiles))

resol_factor = args.resol_factor
print("resol_factor= "+ str(resol_factor))

Nside = args.Nside
print("Nside= "+ str(Nside))

# Access the parsed values
rangeAng = args.rangeAng
print("Parsed values=", range)

downSample = args.downSample
print("Down Sample= "+ str(downSample))



# # parameters

# folderpath = args.folderpath
# # folderpath =  '/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# folderpath =  "/home/asus/Dropbox/extras/storage/graham/ht/data_cps512_512/4Mpc_2048c_1024p_zi63_nowakem/"

# # folderpath =  '/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/'
# # folderpath =  '/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_nowakem/'
# # folderpath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/half_lin_cutoff_half_tot_pert_nvpw/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/"
# # folderpath =  "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_384/14--11-11pv_1.5708-0.84153-3.0434ra/2dproj/dm/"
# # folderpath =  '/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE8/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/data/1lf_1rf/NSIDE_8/anglid_1/-4-11--24pv_0.10211--0.62099-0.7854ra/2dproj/dm/'
# # folderpath =  '/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE4_tst/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/sample1001/half_lin_cutoff_half_tot_pert_nvpw/data/1lf_1rf/NSIDE_4/anglid_86/-15--14--17pv_1.4033-0.79983-5.1051ra/2dproj/dm/'
# # folderpath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/'
# # folderpath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5029/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_21/-232--108-113pv_0.62237--1.5029-4.4506ra/2dproj/dm/'
# # folderpath =  '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_86/-153--150--178pv_1.4033-0.79983-5.1051ra/2dproj/dm/'
# print("File Path in = "+ str(folderpath))

# sample = args.sample
# sample =  'sample5004'
# # sample =  'sample5001'
# # sample =  'sample5029'
# print("sample = "+ str(sample))

# sampleC = args.sampleC
# # sampleC =  '/half_lin_cutoff_half_tot_pert_nvpw_v0p6'
# print("sampleC = "+ str(sampleC))

# path_out = args.path_out
# # path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"
# path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps512_512_hpx_2d_NSIDE4_stat/4Mpc_2048c_1024p_zi63_nowakem/"

# # path_out = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat_tst/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"
# # path_out = "/home/asus/Dropbox/extras/storage/graham/small_res/data_cps12_48_hpx_2d_NSIDE4_tst/plots/64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m/"
# # path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512/plots_4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# # path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4_tst/plots_1/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# # path_out = '/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4/plots_3/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/'
# print("Path out = "+ str(path_out))

# redshift = args.redshift
# # redshift = '63'
# # redshift = '5'
# # redshift = '10'
# redshift = '3'
# print("redshift= "+ str(redshift))



# Nmesh =  ast.literal_eval(args.Nmesh)
# # Nmesh = [48,48,48]
# # Nmesh = [48,48,12]
# # Nmesh = [512,512,32]
# Nmesh = [512,512,512]
# print("Nmesh= "+ str(Nmesh))


# BoxSize = args.BoxSize
# print("BoxSize= "+ str(BoxSize))

# nfiles = args.nfiles
# print("nfiles= "+ str(nfiles))

# resol_factor = args.resol_factor
# # resol_factor = 1
# # resol_factor = 2
# print("resol_factor= "+ str(resol_factor))

# Nside = args.Nside
# Nside = 4
# print("Nside= "+ str(Nside))


# # Access the parsed values
# rangeAng = args.rangeAng
# # rangeAng = parse_range('86')
# # rangeAng = parse_range('21')
# # rangeAng = parse_range('1')
# # rangeAng = parse_range('2')
# # rangeAng = parse_range('1-2,86')
# rangeAng = parse_range('')
# print("Parsed values=", rangeAng)


# downSample = args.downSample
# downSample = 4
# print("Down Sample= "+ str(downSample))

#%%

import numpy as np
import fcntl
import os


# def write_array_at_line(filepath, line_number, array_data):
#     # Convert the array data to a string with tab-separated values
#     line_content = '\t'.join(array_data.astype(int).astype(str)) + '\n'
    
#     # Read existing lines or initialize empty if file does not exist
#     try:
#         with open(filepath, 'r') as file:
#             lines = file.readlines()
#     except FileNotFoundError:
#         lines = []
    
#     # Ensure the list has enough lines by adding blank lines as needed
#     while len(lines) < line_number:
#         lines.append('\n')
    
#     # Update the specific line with the array content
#     lines[line_number - 1] = line_content
    
#     # Write the updated content back to the file
#     with open(filepath, 'w') as file:
#         file.writelines(lines)

def write_array_at_line(filepath, line_number, array_data):
    # If line_number is None (remember this is angleId), write to the first line
    if line_number is None:
        line_number = 1
        
    # Convert array data to integers, then to tab-separated string format
    line_content = '\t'.join(array_data.astype(int).astype(str)) + '\n'
    
    # Read the existing lines or initialize an empty list if the file doesn't exist
    try:
        with open(filepath, 'r+') as file:
            # Lock the file for writing
            fcntl.flock(file, fcntl.LOCK_EX)

            lines = file.readlines()

            # Ensure there are enough lines by appending blank lines if necessary
            while len(lines) < line_number:
                lines.append('\n')

            # Replace the specific line with the array content
            lines[line_number - 1] = line_content

            # Move the file pointer to the beginning and overwrite file with new content
            file.seek(0)
            file.writelines(lines)
            file.truncate()

            # Release the file lock
            fcntl.flock(file, fcntl.LOCK_UN)

    except FileNotFoundError:
        with open(filepath, 'w') as file:
            fcntl.flock(file, fcntl.LOCK_EX)
            
            # Write blank lines up to the required line, then add array content
            lines = ['\n'] * (line_number - 1)
            lines.append(line_content)
            file.writelines(lines)
            fcntl.flock(file, fcntl.LOCK_UN)
        
# Ensure the directory exists
os.makedirs(path_out, exist_ok=True)        

def downsample_block_sum(mesh, l: int):
    '''

    Parameters
    ----------
    mesh : TYPE
        DESCRIPTION.
    l : int
        the donwsample cell number to sum up. Since we are working with multiple of 2, it should also be of this kind

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    mesh = np.asarray(mesh)
    Nx, Ny, Nz = mesh.shape

    assert Nx % l == 0 and Ny % l == 0 and Nz % l == 0, \
        "Each dimension must be divisible by l (or use the pad/crop version)."

    Mx, My, Mz = Nx // l, Ny // l, Nz // l
    
    # out = mesh.reshape(Mx, l, My, l, Mz, l).sum(axis=(1, 3, 5))

    return mesh.reshape(Mx, l, My, l, Mz, l).sum(axis=(1, 3, 5)) / (l**3)   # average; remove this line if you want sum




#%%

import sys

import numpy as np



path_analy = os.getcwd() +'/' 
sys.path.append(path_analy+'read')
import Read_slices



        
for AngId in (rangeAng if rangeAng else [None]):
# for AngId in range(1,2):    
    if AngId is None:
        print("Looking at the only Angle")
        folder_path = folderpath + sample + sampleC + '/data/' 
        subfolder_name = os.listdir(folder_path)[0]
    else:
        print("Looking at angle id = ", AngId)
        # mesh = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)
        folder_path = folderpath + sample + sampleC + '/data/' + "1lf_"+ str(resol_factor) + 'rf/NSIDE_' + str(Nside) + '/anglid_' + str(AngId) + '/'
        
    subfolder_name = os.listdir(folder_path)[0]
    filepath = os.path.join(folder_path, subfolder_name)+'/2dproj/dm/'
    mesh = Read_slices.read_slices_bin(Nmesh,BoxSize,nfiles,filepath,redshift)
    
    mesh = downsample_block_sum(mesh,downSample)
    
    numVoid_array = np.zeros(int(Nmesh[2]/downSample+1))
    sum_void = 0
    
    for slic in range(0,int(Nmesh[2]/downSample)):
    # for slic in range(1,1+1):    
        slice_data = mesh[:,:,slic]
        valuesVoid = slice_data.flatten()
        valuesVoid = valuesVoid[valuesVoid == -1.]
        numVoid = np.count_nonzero(valuesVoid)
        numVoid_array[slic]=numVoid
        sum_void = sum_void + numVoid
    
    numVoid_array[int(Nmesh[2]/downSample)]=sum_void
    
    write_array_at_line(path_out+sample+'_2d_void_z'+redshift+"_down"+str(downSample)+"_stat.txt", AngId, numVoid_array)


    
    
    

#%%


# import numpy as np
# import fcntl

# def write_array_at_line(filepath, line_number, array_data):
#     # Convert array data to integers, then to tab-separated string format
#     line_content = '\t'.join(array_data.astype(int).astype(str)) + '\n'
    
#     # Read the existing lines or initialize an empty list if the file doesn't exist
#     try:
#         with open(filepath, 'r+') as file:
#             # Lock the file for writing
#             fcntl.flock(file, fcntl.LOCK_EX)

#             lines = file.readlines()

#             # Ensure there are enough lines by appending blank lines if necessary
#             while len(lines) < line_number:
#                 lines.append('\n')

#             # Replace the specific line with the array content
#             lines[line_number - 1] = line_content

#             # Move the file pointer to the beginning and overwrite file with new content
#             file.seek(0)
#             file.writelines(lines)
#             file.truncate()

#             # Release the file lock
#             fcntl.flock(file, fcntl.LOCK_UN)

#     except FileNotFoundError:
#         with open(filepath, 'w') as file:
#             fcntl.flock(file, fcntl.LOCK_EX)
            
#             # Write blank lines up to the required line, then add array content
#             lines = ['\n'] * (line_number - 1)
#             lines.append(line_content)
#             file.writelines(lines)
#             fcntl.flock(file, fcntl.LOCK_UN)



# def write_array_at_line(filepath, line_number, array_data):
#     # Convert the array data to a string with tab-separated values
#     line_content = '\t'.join(array_data.astype(int).astype(str)) + '\n'
    
#     # Read existing lines or initialize empty if file does not exist
#     try:
#         with open(filepath, 'r') as file:
#             lines = file.readlines()
#     except FileNotFoundError:
#         lines = []
    
#     # Ensure the list has enough lines by adding blank lines as needed
#     while len(lines) < line_number:
#         lines.append('\n')
    
#     # Update the specific line with the array content
#     lines[line_number - 1] = line_content
    
#     # Write the updated content back to the file
#     with open(filepath, 'w') as file:
#         file.writelines(lines)




# def write_array_at_line(filepath, line_number, array_data):
#     # Convert the array data to a string with tab-separated values
#     line_content = '\t'.join(array_data.astype(int).astype(str)) + '\n'

#     # Read existing lines or create an empty list if the file doesn't exist
#     try:
#         with open(filepath, 'r') as file:
#             lines = file.readlines()
#     except FileNotFoundError:
#         lines = []

#     # Ensure the list has enough lines by adding blank lines as needed
#     while len(lines) < line_number:
#         lines.append('\n')

#     # Update the specific line with the array content
#     lines[line_number - 1] = line_content

#     # Write the updated content back to the file
#     with open(filepath, 'w') as file:
#         file.writelines(lines)
        

        
        
#%%



# file = open(path_out+sample+'_2d_void_z'+redshift+"_stat.txt", "w")

# write_array_at_line(path_out+sample+'_2d_void_z'+redshift+"_stat.txt", AngId, numVoid_array)


#%%


# # explore data

# # import matplotlib
# # save=[]
# # if save is None:
# #     matplotlib.use('Qt5Agg')    # to show figures on desktop
# # else:
# #     matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)

# from matplotlib import pyplot as plt
# import numpy as np


# slice_data = slice_data = mesh[:,:,182]
# plt.figure()    
#     # plt.imshow(mesh.preview(axes=[0,2]))
#     # plt.imshow(np.log10(1+mesh.preview(axes=[0,1])))
# plt.imshow(np.log10(+1+slice_data))










    