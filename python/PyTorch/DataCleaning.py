#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:43:51 2024

@author: asus
"""



# load data 
import os
from PIL import Image
import numpy as np

#just to debug
import matplotlib.pyplot as plt

#root directory
folder_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"

#function returns all figures paths inside folder_path including subfolders
def list_all_files(folder_path):
    all_files = []
    
    # os.walk generates the file names in a directory tree, including all subdirectories
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            # Construct full file path
            file_path = os.path.join(root, file)
            all_files.append(file_path)
    
    return all_files



# Example usage
# folder_path = "/path/to/your/folder"
files_list = list_all_files(folder_path)

# # It is possible to separate in classes afterward
# classes = sorted(set(os.path.basename(os.path.dirname(path)) for path in files_list))

# filename = files_list[0]

# # Open the image file
# with Image.open(filename) as img:
#     # Convert the image to grayscale (black and white)
#     gray_img = img.convert("L")    

#     # Convert the grayscale image to a numpy array
#     intensity_array = np.array(gray_img)    

# # plt.imshow(gray_img, cmap='gray')



# # # debug
# values = intensity_array.flatten()
# values = values[values != 38.]

# # number of outliers
# num_out1 = np.count_nonzero(values <= 60.)

# a = values.min()


# minim_log = np.floor(np.log10(values.min()))
# maxim_log = np.log10(values.max())
# num_per_log10 = 100
# numb_bins = int((maxim_log-minim_log)*num_per_log10)
# bins = np.logspace(minim_log, maxim_log, numb_bins)    
# hist2d = plt.hist(values, bins=bins)
# # plt.xscale('log')
# plt.yscale('log')












# # dictionalry with the number of outliers for each file, for debbuging, since 
# # it will store all values
# outlier_lists = {}

# list with the files that passes the threashold test
filtered_filenames = []
# dic with the files that passes the threashold test (only for visualization)
filtered_filenames_dic = {}

# Loop through all folder_path
for file_path in files_list:
      
    # Open the image file
    with Image.open(file_path) as img:
        # Convert the image to grayscale (black and white)
        gray_img = img.convert("L")
        
        # Convert the grayscale image to a numpy array
        intensity_array = np.array(gray_img)        
        
        # obtain values on the intesnity map
        values = intensity_array.flatten()
        # remove the black ticks and black borders 
        values = values[values != 38.]
        # number of outliers
        num_out = np.count_nonzero(values <= 60.)
        
        if num_out<=4:
            filtered_filenames.append(file_path)
            filtered_filenames_dic[file_path] = num_out
        # # Store the intensity array in the dictionary, for debugging
        # outlier_lists[file_path] = num_out        


# # Split the image into individual channels
# r, g, b = img.split()

# # Plot each channel
# fig, ax = plt.subplots(1, 3, figsize=(12, 4))

# # Display the Red channel
# ax[0].imshow(r, cmap='gray')
# ax[0].set_title('Red Channel')
# ax[0].axis('off')  # Hide axes ticks
# # Convert the grayscale image to a numpy array
# intensity_arrayR = np.array(r)  

# # Display the Green channel
# ax[1].imshow(g, cmap='gray')
# ax[1].set_title('Green Channel')
# ax[1].axis('off')
# intensity_arrayG = np.array(g)  

# # Display the Blue channel
# ax[2].imshow(b, cmap='gray')
# ax[2].set_title('Blue Channel')
# ax[2].axis('off')
# intensity_arrayB = np.array(b)  

# # Show the plot
# plt.show()
