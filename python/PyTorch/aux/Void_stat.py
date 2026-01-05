#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 13:27:17 2024

@author: asus
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 09:56:55 2024

@author: Disrael
"""


import numpy as np

#%%

# colect data

# File path and specifications
path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat/void/"
wake_spec = ["4Mpc_2048c_1024p_zi63_nowakem/","4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"]

rang=range(3001,3010+1)
# rang=range(5001,5100+1)


n_angle = 96
slices = 33
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
    filename_nowake = path + wake_spec[0] + f"sample{simul}_2d_void_z3_stat.txt"
    filename_wake =   path+wake_spec[1]+"sample"+str(simul)+"_2d_void_z3_stat.txt" 
    
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

# plot data


import matplotlib.pyplot as plt

# Example 3D array
# Replace this with your actual 3D array data

num_per_log10 = 10
num_per_1 = 2



# Flatten the 3D array into a 1D array
all_data_nowake_flat = all_data_nowake.flatten()
all_data_wake_flat = all_data_wake.flatten()


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
# plt.xscale('log')
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
# plt.xscale('log')
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of 3D Array Values wake')
plt.show()
        
        