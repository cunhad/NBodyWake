#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 09:56:55 2024

@author: Disrael
"""


import numpy as np



# File path and specifications
path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_stat_2dc1l1_3dc1l1/"
wake_spec = ["4Mpc_2048c_1024p_zi63_nowakem/","4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"]

rang=range(3001,3010+1)
# rang=range(5001,5100+1)


n_angle = 96
slices = 32
# slices_dp = 8

# Initialize an empty list to store the flattened data
all_data_nowake = np.empty((len(rang), n_angle, slices))
all_data_wake = np.empty((len(rang), n_angle, slices))



# for simul in range(3001,3010):
# for simul in range(5001,5100):
# for simul in range(3001,3001+1):
for i, simul in enumerate(rang):    
    # print("sample"+str(simul))
    # filename_nowake = path+wake_spec[0]+"sample"+str(simul)+"_2ds4t3_curv_z3_stat.txt" 
    filename_nowake = path + wake_spec[0] + f"sample{simul}_2ds4t3_curv_z3_stat.txt"
    filename_wake =   path+wake_spec[1]+"sample"+str(simul)+"_2ds4t3_curv_z3_stat.txt" 
    
    try:
        data_array_nowake = np.loadtxt(filename_nowake, delimiter='\t')
        all_data_nowake[i, :, :] = data_array_nowake
        # Check if the array has the expected dimensions
        if data_array_nowake.shape != (n_angle, slices):
            print(f"The array no wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_nowake.shape}.")
            break
        # Flatten the data array to 1D and append to the list
        # all_data_nowake.extend(data_array_nowake.flatten())
        
    except Exception as e:
        print(f"Error loading file {filename_nowake}: {e}")
        break
    
    try:
        data_array_wake = np.loadtxt(filename_wake, delimiter='\t')
        all_data_wake[i, :, :] = data_array_wake
        # Check if the array has the expected dimensions
        if data_array_wake.shape != (n_angle, slices):
            print(f"The array wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_wake.shape}.")
            break
        # Flatten the data array to 1D and append to the list
        # all_data_wake.extend(data_array_wake.flatten())
        
    except Exception as e:
        print(f"Error loading file {filename_nowake}: {e}")
        break
    

# # Convert the list to a numpy array if needed
# all_data_array_nowake = np.array(all_data_nowake)
# all_data_array_wake = np.array(all_data_wake)




# #%%

# # Simplified test for a single file within a loop
# for simul in [3006]:  # Only test with one file
#     filename_nowake = path + wake_spec[0] + "sample" + str(simul) + "_2ds4t3_curv_z3_stat.txt"
#     try:
#         data_array_nowake = np.loadtxt(filename_nowake, delimiter='\t')
#         print("File loaded successfully with shape:", data_array_nowake.shape)
#     except Exception as e:
#         print(f"Error loading file {filename_nowake}: {e}")
        
        
#%%



import matplotlib.pyplot as plt

# Example 3D array
# Replace this with your actual 3D array data

num_per_log10 = 10
num_per_1 = 2


# Difference data

all_data_difference = all_data_wake - all_data_nowake

# Flatten the 3D array into a 1D array
all_data_nowake_flat = all_data_nowake.flatten()
all_data_wake_flat = all_data_wake.flatten()
all_data_difference_flat = all_data_difference.flatten()


# find a commom x-range:
minim = min(all_data_nowake_flat.min(),all_data_wake_flat.min())
maxim = max(all_data_nowake_flat.max(), all_data_wake_flat.max())
numb_bins = int((maxim-minim)*num_per_1)
bins = np.linspace(minim, maxim, numb_bins)   

# minim_log = np.log10(all_data_nowake_flat.min())
# maxim_log = np.log10(all_data_nowake_flat.max())
# numb_bins = int((maxim_log-minim_log)*num_per_log10)
# bins = np.logspace(minim_log, maxim_log, numb_bins)   

# minim = all_data_nowake_flat.min()
# maxim = all_data_nowake_flat.max()
# numb_bins = int((maxim-minim)*num_per_1)
# bins = np.linspace(minim, maxim, numb_bins)   

# Plot the histogram
# hist3d = plt.hist(values, bins=bins)
plt.figure()  
hist3d = plt.hist(all_data_nowake_flat, bins=bins, color='skyblue', edgecolor='black')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of 3D Array Values no wake')
plt.show()
        
        

# Flatten the 3D array into a 1D array

# minim_log = np.log10(all_data_nowake_flat.min())
# maxim_log = np.log10(all_data_nowake_flat.max())
# numb_bins = int((maxim_log-minim_log)*num_per_log10)
# bins = np.logspace(minim_log, maxim_log, numb_bins)   

# minim = all_data_wake_flat.min()
# maxim = all_data_wake_flat.max()
# numb_bins = int((maxim-minim)*num_per_1)
# bins = np.linspace(minim, maxim, numb_bins)   

# Plot the histogram
# hist3d = plt.hist(values, bins=bins)
plt.figure()  
hist3d = plt.hist(all_data_wake_flat, bins=bins, color='skyblue', edgecolor='black')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of 3D Array Values wake')
plt.show()
        

#%%

# Define x-range:
# minim = min(all_data_difference_flat)
# maxim = max(all_data_difference_flat)
minim = 0
maxim = 100
num_per_log10 = 10
num_per_1 = 2
numb_bins = int((maxim-minim)*num_per_1)
bins = np.linspace(minim, maxim, numb_bins)

plt.figure()  
hist3d = plt.hist(all_data_difference_flat, bins=bins, color='skyblue', edgecolor='black')
# plt.xscale('log')
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of 3D Array Values wake no wake Diff')
plt.show()


#%%


# max_value = np.max(all_data_wake)
# max_index_flat = np.argmax(all_data_wake)
# max_indices = np.unravel_index(max_index_flat, all_data_wake.shape)

# flat = all_data_wake.ravel()
# flat = flat[~np.isnan(flat)]

# top2 = np.partition(flat, -2)[-2:]
# second_max = top2.min()

# print the k highest
k = 200
# threshold = 3
threshold = np.inf

flat = all_data_wake.ravel()
# flat = all_data_nowake.ravel()
# flat = all_data_difference.ravel()
# Mask invalid values
mask = (~np.isnan(flat)) & (flat < threshold)
if np.count_nonzero(mask) < k:
        raise ValueError("Not enough values below threshold to extract top-k.")
flat_masked = np.where(mask, flat, -np.inf)

# # Handle NaNs
# flat_clean = np.where(np.isnan(flat), -np.inf, flat)

idxs = np.argpartition(flat_masked, -k)[-k:]
idxs = idxs[np.argsort(flat_masked[idxs])[::-1]]

values = flat_masked[idxs]
indices = [np.unravel_index(i, all_data_wake.shape) for i in idxs]

print(values)
print(indices)
        