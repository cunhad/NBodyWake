#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 09:42:19 2024

@author: asus
"""

import os


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



# Load datasets

# folder_path = "/scratch/cunhad/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"
# folder_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50_test2/"
folder_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50_test2/"
files_list = list_all_files(folder_path)



from PIL import Image
import numpy as np

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

import random


def select_valid_classes(classes_dict,threshold):

    keys = list(classes_dict.keys())
    
 
    
    indices = list(range(len(classes_dict)))
    random.shuffle(indices)
    
    selected_samples = []
    current_sum = 0
    selected_indices = []
    selected_samples_complement = keys

    

    
    for index in indices:
        # print(index)
        # print(current_sum > threshold)
        if current_sum > threshold:
            # print(index)
            break
        selected_samples.append(keys[index])
        # b=classes_dict[keys[index]['count']
        current_sum += classes_dict[keys[index]]['count']
        selected_indices.append(index)
        
    selected_indices.sort(reverse=True)
    # Remove elements at the specified indices
    for index in selected_indices:
        del selected_samples_complement[index]
        
    return selected_samples,current_sum,selected_samples_complement


def file_list(folder_path,VALID_RATIO):
    files_list = list_all_files(folder_path)
    filtered_filenames,filtered_filenames_dic = number_of_outliers(files_list)
    classes_dict = classify_files(filtered_filenames)

    # Define the proportion or number of items in each set

    n_train_examples = int(len(filtered_filenames) * VALID_RATIO)
    n_valid_examples = len(filtered_filenames) - n_train_examples

    selected_indices,current_sum,selected_indices_complement=select_valid_classes(classes_dict,n_valid_examples)
    
    # Get all associated elements in the dictionary for the given keys list
    val_files = [classes_dict[key]["filenames"] for key in selected_indices if key in classes_dict]
    # Flatten the list of lists into a 1D list using list comprehension
    val_files = [item for sublist in val_files for item in sublist]


    # Get the complementary elements (keys not in the list)
    # test_files = {key: sample_dict[key] for key in sample_dict if key not in keys_list}
    test_files =  [classes_dict[key]["filenames"] for key in classes_dict if key not in selected_indices]
    test_files = [item for sublist in test_files for item in sublist]

    
    return val_files,test_files,selected_indices,selected_indices_complement


VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation

valid_data_, test_train_data_,selected_indices_val,selected_indices_test_train = file_list(folder_path,VALID_RATIO)


#%%
def count_files_by_subdirectory(file_list, subdirectories):
    counts = {subdir: 0 for subdir in subdirectories}
    
    for file_path in file_list:
        for subdir in subdirectories:
            if subdir in file_path:
                counts[subdir] += 1
    
    return counts

# count the excess
subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in valid_data_])))
counts = count_files_by_subdirectory(valid_data_, subdirectories)
print("Valid")
print(counts)

# count the excess
subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in test_train_data_])))
counts = count_files_by_subdirectory(test_train_data_, subdirectories)
print("Test and Train")
print(counts)



#%%




def isClean(filename):
    # Open the image file
    with Image.open(filename) as img:
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
            return True
        else:
            return False     


import re


def balanced_list_of_files(files_list_original,files_list_complement,files_list_to_balance,selected_indices):
    files_list_only_names = [file.split("/")[-1] for file in files_list_original]
    
    # count the excess
    subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_to_balance])))
    counts = count_files_by_subdirectory(files_list_to_balance, subdirectories)
    wake_excess = counts[subdirectories[1]]-counts[subdirectories[0]]
    
    pattern = re.compile(r'sample(\d+)-anglid_\d+-2dproj_z3_ts\d+_sl\d+\.png')

    
    # If there are more wakes, add no wakes
    if wake_excess> 0:
    
        # Obtain the elements of the complementar list
        files_list_balance_ = [file for file in files_list_complement if subdirectories[1] not in file]
        files_list_balance__ = [file for file in files_list_balance_ if file.split("/")[-1] not in files_list_only_names]
        files_list_balance___ = [element for element in files_list_balance__ if pattern.search(element).group(1) in selected_indices]
        
        # For each missing smaple, test if pass the noise test. if yes, add it
        files_list_balance____ = []
        for i in range(len(files_list_balance___)):
            
            # Choose a random element from the list
            random_element = random.choice(files_list_balance___)
            
            # Remove the chosen element from the list
            files_list_balance___.remove(random_element)
            if isClean(random_element):
                files_list_balance____.append(random_element)
                wake_excess -= 1
                if wake_excess == 0:
                    break

        files_list_balanced = files_list_to_balance + files_list_balance____
        
    else:
        for i in range(abs(wake_excess)):
            
            # Obtain the list of no wake files
            files_list_no_wake = [file for file in files_list_to_balance if subdirectories[1] not in file]
            # Choose a random element from the list
            random_element = random.choice(files_list_no_wake)
            # Remove the chosen element from the list
            files_list_to_balance.remove(random_element)
            
        
        files_list_balanced = files_list_to_balance + files_list_to_balance
        
    return files_list_balanced
            



# folder_path_balance = "/scratch/cunhad/data_cps32_512_hpx_2d_NSIDE4_figs_thr40/"
folder_path_balance = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"
files_list_balance  = list_all_files(folder_path_balance)


files_list_balanced_val = balanced_list_of_files(files_list,files_list_balance,valid_data_,selected_indices_val)
files_list_balanced_test_train = balanced_list_of_files(files_list,files_list_balance,test_train_data_,selected_indices_test_train)

#%%


# import re


# # folder_path_balance = "/scratch/cunhad/data_cps32_512_hpx_2d_NSIDE4_figs_thr40/"
# folder_path_balance = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"
# files_list_balance  = list_all_files(folder_path_balance)



# files_list_original = files_list
# files_list_complement = files_list_balance
# files_list_to_balance = test_train_data_
# selected_indices = selected_indices_test_train

# files_list_only_names = [file.split("/")[-1] for file in files_list_original]

# # count the excess
# subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_to_balance])))
# counts = count_files_by_subdirectory(files_list_to_balance, subdirectories)
# wake_excess = counts[subdirectories[1]]-counts[subdirectories[0]]


# pattern = re.compile(r'sample(\d+)-anglid_\d+-2dproj_z3_ts\d+_sl\d+\.png')


# # If there are more wakes, add no wakes
# if wake_excess> 0:
#     # Obtain the elements of the complementar list
#     files_list_balance_ = [file for file in files_list_complement if subdirectories[1] not in file]
#     files_list_balance__ = [file for file in files_list_balance_ if file.split("/")[-1] not in files_list_only_names]
#     files_list_balance___ = [element for element in files_list_balance__ if pattern.search(element).group(1) in selected_indices]


#     # For each missing smaple, test if pass the noise test. if yes, add it
#     files_list_balance____ = []
#     for i in range(len(files_list_balance___)):
#         print(wake_excess, " try")
#         # Choose a random element from the list
#         random_element = random.choice(files_list_balance___)
        
#         # Remove the chosen element from the list
#         files_list_balance___.remove(random_element)
#         if isClean(random_element):
#             files_list_balance____.append(random_element)
#             wake_excess -= 1
#             if wake_excess == 0:
#                 break
#             print( "add")
#         else:
#             print("no add")
    
#     files_list_balanced = files_list_to_balance + files_list_balance____



# # count the excess
# subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_balanced])))
# counts = count_files_by_subdirectory(files_list_balanced, subdirectories)
# print("Test and Train")
# print(counts)


# # count the excess
# subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_balance____])))
# counts = count_files_by_subdirectory(files_list_balance____, subdirectories)
# print("Test and Train")
# print(counts)


#%%

# count the excess
subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_balanced_val])))
counts = count_files_by_subdirectory(files_list_balanced_val, subdirectories)
print("Valid")
print(counts)

# count the excess
subdirectories = sorted(list(set([os.path.basename(os.path.dirname(file_path)) for file_path in files_list_balanced_test_train])))
counts = count_files_by_subdirectory(files_list_balanced_test_train, subdirectories)
print("Test and Train")
print(counts)



