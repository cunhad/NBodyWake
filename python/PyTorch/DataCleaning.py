#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:43:51 2024

@author: asus
"""

VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation

#root directory
folder_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"



# load data 
import os
from PIL import Image
import numpy as np


import random


# #just to debug
# import matplotlib.pyplot as plt

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





#%%




# Function that returns the number of outliers (empty cell)

def number_of_outliers(files_list):

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
        
    return filtered_filenames,filtered_filenames_dic


#%%




def classify_files(file_list):
    import re
    from collections import defaultdict

    # Regular expression to capture the class information from the file path
    pattern = re.compile(r'/sample(\d+)-')

    # Dictionary to store class information
    classes = defaultdict(lambda: {"count": 0, "filenames": []})

    for file in file_list:
        match = pattern.search(file)
        if match:
            class_num = match.group(1)
            classes[class_num]["count"] += 1
            classes[class_num]["filenames"].append(file)
    
    return classes





#%%



def select_valid_classes(classes_dict,threshold):

    keys = list(classes_dict.keys())
    
 
    
    indices = list(range(len(classes_dict)))
    random.shuffle(indices)
    
    selected_indices = []
    current_sum = 0
    

    
    for index in indices:
        # print(index)
        # print(current_sum > threshold)
        if current_sum > threshold:
            # print(index)
            break
        selected_indices.append(keys[index])
        # b=classes_dict[keys[index]['count']
        current_sum += classes_dict[keys[index]]['count']
        
    return selected_indices,current_sum





# val_files = classes_dict[selected_indices]['count']


def file_list(folder_path,VALID_RATIO):
    files_list = list_all_files(folder_path)
    filtered_filenames,filtered_filenames_dic = number_of_outliers(files_list)
    classes_dict = classify_files(filtered_filenames)

    # Define the proportion or number of items in each set

    n_train_examples = int(len(filtered_filenames) * VALID_RATIO)
    n_valid_examples = len(filtered_filenames) - n_train_examples

    selected_indices,current_sum=select_valid_classes(classes_dict,n_valid_examples)
    
    # Get all associated elements in the dictionary for the given keys list
    val_files = [classes_dict[key]["filenames"] for key in selected_indices if key in classes_dict]
    # Flatten the list of lists into a 1D list using list comprehension
    val_files = [item for sublist in val_files for item in sublist]


    # Get the complementary elements (keys not in the list)
    # test_files = {key: sample_dict[key] for key in sample_dict if key not in keys_list}
    test_files =  [classes_dict[key]["filenames"] for key in classes_dict if key not in selected_indices]
    test_files = [item for sublist in test_files for item in sublist]

    
    return val_files,test_files