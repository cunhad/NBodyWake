#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:41:50 2024

@author: Disrael
"""

#%%


# parser
# Run example:
# python WakeDetection_Pytorch_CNN_TranfLearn_term --batch_size=32 --num_workers=0

def parse_range(value):
    numbers = set()
    parts = value.split(',')

    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            numbers.update(range(start, end + 1))
        else:
            numbers.add(int(part))

    return sorted(numbers)

# def parse_range(value):
#     numbers = set()
#     parts = value.split(',')

#     for part in parts:
#         if '-' in part:
#             start, end = map(int, part.split('-'))
#             numbers.update(range(start, end + 1))
#         else:
#             numbers.add(int(part))

#     sorted_numbers = sorted(numbers)

#     # Convert sorted numbers into ranges
#     ranges = []
#     start = sorted_numbers[0]
#     prev = start

#     for num in sorted_numbers[1:]:
#         if num != prev + 1:  # If gap detected, close current range
#             ranges.append(range(start, prev + 1))
#             start = num
#         prev = num

#     # Append the last range
#     ranges.append(range(start, prev + 1))

#     return ranges



import argparse

parser = argparse.ArgumentParser(description='efficientnet_b7 wake classification model')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--path_data', type=str, help='')
parser.add_argument('--path_void', type=str, help='')
parser.add_argument('--path_WakeSignal', type=str, help='')

parser.add_argument('--n_angle', type=int, default=96, help='')
parser.add_argument('--rangeSampl', type=parse_range, default='5001-5100', help="Range of samples ids (e.g., '5001-5100,5200,5211')")

parser.add_argument('--slices_void', type=int, default=33, help='')
parser.add_argument('--slices_signal', type=int, default=32, help='')
parser.add_argument('--percentage_positiveWakeSig', type=float, default=10, help='')
parser.add_argument('--validation_fraction', type=float, default=0.1, help='')
parser.add_argument('--train_tt_fraction', type=float, default=0.8, help='')
parser.add_argument('--wake_top_percentage', type=float, default=20, help='')
parser.add_argument('--void_percentage', type=float, default=0, help='')


parser.add_argument('--batch_size', type=int, default=32, help='')
parser.add_argument('--num_workers', type=int, default=0, help='')
parser.add_argument('--num_epochs', type=int, default=10, help='')



args = parser.parse_args()

# parameters

path_data = args.path_data
path_data = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs/"
print("File Path in = "+ str(path_data))

path_void = args.path_void
path_void = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat/void/"
print("File Path Void in = "+ str(path_void))

path_WakeSignal = args.path_WakeSignal
path_WakeSignal = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_stat_2dc1l1_3dc1l1/"
print("File Path Wake Signal in = "+ str(path_WakeSignal))




n_angle = args.n_angle
# num_epochs = 1
print("angles each sample = "+ str(n_angle))

# Access the parsed values
rang = args.rangeSampl
# rangeAng = parse_range('5001-5101')
print("Samples = ", rang)




slices_void = args.slices_void
# slices_signal = 33
print("Slices for void = "+ str(slices_void))

slices_signal = args.slices_signal
# slices_signal = 32
print("Slices for wake signal = ", slices_signal)

percentage_positiveWakeSig = args.percentage_positiveWakeSig
# percentage_positiveWakeSig = 10
print("Percentage positive Wake Signal = ", percentage_positiveWakeSig)

validation_fraction = args.validation_fraction
# validation_fraction = 0.1
print("Validation fraction of total data = ", validation_fraction)

train_tt_fraction = args.train_tt_fraction
# train_tt_fraction = 0.1
print("Train fraction of train + test data = ", train_tt_fraction)


wake_top_percentage = args.wake_top_percentage
# wake_top_percentage = 10
print("Percentage top signal wake = ", wake_top_percentage)

void_percentage = args.void_percentage
# void_percentage = 0.1
print("Percentage lower voids wake = ", void_percentage)




batch_size =  args.batch_size
batch_size =  4
print("Batch size = "+ str(batch_size))

num_workers = args.num_workers
# num_workers = 0
print("Num of CPU workers = "+ str(num_workers))

num_epochs = args.num_epochs
# num_epochs = 1
print("Num epochs = "+ str(num_epochs))


# # batch_size = 32
# TRAIN_RATIO = 0.8       #fraction of total dataset that will go to train+validation
# VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation
OUTPUT_DIM = 1          # 2 classes for classification labels
# SEED = 1234
pretrained_size = 512


wake_spec = ["4Mpc_2048c_1024p_zi63_nowakem/","4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"]



#%%

# Import the necessary libraries:

# Torch stuff
    
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset,random_split
from torchvision import transforms, datasets
import torch.utils.data as data

# from efficientnet_pytorch import EfficientNet
import torchvision.models as models

# Training function def
import copy 

# data load

import os
from PIL import Image

# set the random seeds

import random

# math

import numpy as np


# Import the spliter function from DataCleaning.py'

# from DataCleaning import file_list,list_all_files,balanced_list_of_files

from Data_preSelection  import  *

# keep track of true positives, false positives, true negatives, and false negatives for each class
from sklearn.metrics import confusion_matrix


# To count the number of images in each label
from collections import Counter





#%%

# Define the data transformations:

    
pretrained_means = [0.485, 0.456, 0.406]
pretrained_stds = [0.229, 0.224, 0.225]

train_transforms = transforms.Compose([
                           transforms.Resize(pretrained_size),
                           transforms.RandomRotation(5),
                           transforms.RandomHorizontalFlip(0.5),
                           transforms.RandomCrop(pretrained_size, padding=10),
                           transforms.ToTensor(),
                           transforms.Normalize(mean=pretrained_means,
                                                std=pretrained_stds)
                       ])

test_transforms = transforms.Compose([
                           transforms.Resize(pretrained_size),
                           transforms.ToTensor(),
                           transforms.Normalize(mean=pretrained_means,
                                                std=pretrained_stds)
                       ])





#%%

# Define the Custom Dataset Class

class CustomImageDataset(Dataset):
    def __init__(self, file_list, transform=None):
        """
        Args:
            file_list (list): List of file paths to be included in the dataset.
            transform (callable, optional): Optional transform to be applied on a sample.
        """
        self.file_list = file_list
        self.transform = transform

        # Extract classes from the parent of the immediate parent directory
        self.classes = list(set([os.path.basename(os.path.dirname(os.path.dirname(file_path))) for file_path in file_list]))
        self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}
        self.labels = [self.class_to_idx[os.path.basename(os.path.dirname(os.path.dirname(file_path)))] for file_path in file_list]

    def __len__(self):
        return len(self.file_list)

    def __getitem__(self, idx):
        img_path = self.file_list[idx]
        image = Image.open(img_path).convert("RGB")

        # Extract class name from the parent of the immediate parent directory
        class_name = os.path.basename(os.path.dirname(os.path.dirname(img_path)))
        label = self.class_to_idx[class_name]

        if self.transform:
            image = self.transform(image)

        return image, label


#%%


# Load datasets


all_data_nowake_void, all_data_wake_void = data_void_out(rang, n_angle, slices_void, wake_spec, path_void)
all_data_nowake_signal, all_data_wake_signal = data_signal_out(rang, n_angle, slices_signal, wake_spec, path_WakeSignal)

files_list_all = list_all_files(path_data)
samples_all, anglids_all, tilesizes_all, sliceids_all, wake_infos_all = extract_info(files_list_all)

range_start = rang[0]
data_void = extract_stat(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_void, all_data_wake_void, range_start)
data_signal = extract_stat(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, range_start)
data_signal_diff = extract_stat_diff(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, range_start)

# split val and trainTest

files_list_validation, files_list_trainTest = split_unique_samples(samples_all, files_list_all, validation_fraction)

samples_val, anglids_val, tilesizes_val, sliceids_val, wake_infos_val = extract_info(files_list_validation)
range_start = rang[0]
data_void_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_void, all_data_wake_void, range_start)
data_signal_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, range_start)
data_signal_diff_val = extract_stat_diff(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, range_start)

samples_tt, anglids_tt, tilesizes_tt, sliceids_tt, wake_infos_tt = extract_info(files_list_trainTest)
range_start = rang[0]
data_void_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_void, all_data_wake_void, range_start)
data_signal_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, range_start)
data_signal_diff_tt = extract_stat_diff(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, range_start)



# # Select top wake values, keep void_percentage = 50%
# selected_files, selected_positions, selected_signal_diff = select_extreme_files(
#     data_signal_diff, wake_infos_all, files_list_all, data_void,
#     wake_top_percentage, void_percentage
# )




# validation

# Select top 50% wake values, keep void_percentage = 50%
selected_files_val, selected_positions_val, selected_signal_diff_val = select_extreme_files(
    data_signal_diff_val, wake_infos_val, files_list_validation, data_void_val,
    wake_top_percentage, void_percentage
)

find_extreme_files(selected_signal_diff_val, selected_files_val)



# train and test

# Select top 50% wake values, keep void_percentage = 50%
selected_files_tt, selected_positions_tt, selected_signal_diff_tt = select_extreme_files(
    data_signal_diff_tt, wake_infos_tt, files_list_trainTest, data_void_tt,
    wake_top_percentage, void_percentage
)


find_extreme_files(selected_signal_diff_tt, selected_files_tt)



valid_data__ = CustomImageDataset(file_list=selected_files_val, transform=test_transforms)
test_train_data__ = CustomImageDataset(file_list=selected_files_tt, transform=test_transforms)



n_train_examples = int(len(test_train_data__) * train_tt_fraction)
n_test_examples = len(test_train_data__) - n_train_examples

train_data_, test_data_ = data.random_split(test_train_data__,
                                            [n_train_examples, n_test_examples])



#%%




# Create data loaders.

valid_dataloader = DataLoader(valid_data__, batch_size=batch_size, shuffle=False, num_workers=num_workers)
test_dataloader = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
train_dataloader = DataLoader(train_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)

for X, y in train_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break
    
for X, y in test_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break  

for X, y in valid_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break  



print(test_train_data__.class_to_idx)   # this is the label dictionary









#%%


# Display the number of images in each class for train and test datasets
def count_class_images(subset, original_dataset):
    labels = [original_dataset.labels[i] for i in subset.indices]
    class_counts = Counter(labels)
    class_names = original_dataset.classes
    for class_idx, count in class_counts.items():
        class_name = class_names[class_idx]
        print(f"Class: {class_name}, Number of images: {count}")

print("\nTraining Data:")
count_class_images(train_data_, test_train_data__)

print("\nTesting Data:")
count_class_images(test_data_, test_train_data__)



# Count the number of images in each class
class_counts = Counter(valid_data__.labels)
# Get the class names
class_names = valid_data__.classes

# Display the number of images in each class
print("\nValidation Data:")
for class_idx, count in class_counts.items():
    class_name = class_names[class_idx]
    print(f"Class: {class_name}, Number of images: {count}")


#%%

# Load the EfficientNetB7 model:

model = models.efficientnet_b7(pretrained=True)

# Freeze all layers
for param in model.parameters():
    param.requires_grad = False
    
    
# Modify the final layer for binary classification
IN_FEATURES = model.classifier[-1].in_features
model.classifier[-1] = nn.Linear(IN_FEATURES, OUTPUT_DIM)   

# Unfreeze the last layer
for param in model.classifier[-1].parameters():
    param.requires_grad = True 
    
# Unfreeze the last-to-last layer
for param in model.classifier[-2].parameters():
    param.requires_grad = True    
    
# Define the loss function and optimizer
criterion = nn.BCEWithLogitsLoss()  # Binary Cross Entropy with Logits Loss
optimizer = optim.Adam(
    [{'params': model.classifier[-1].parameters()}, 
     {'params': model.classifier[-2].parameters()}], 
    lr=0.001
)

    
#%%

# # Training and validation functions
# def train_model(model, criterion, optimizer, num_epochs=25):
#     device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#     model = model.to(device)

#     best_model_wts = copy.deepcopy(model.state_dict())
#     best_acc = 0.0

#     for epoch in range(num_epochs):
#         print(f'Epoch {epoch}/{num_epochs - 1}')
#         print('-' * 10)

#         for phase in ['train', 'val']:
#             if phase == 'train':
#                 model.train()  # Set model to training mode
#             else:
#                 model.eval()   # Set model to evaluate mode

#             running_loss = 0.0
#             running_corrects = 0

#             dataloader = train_dataloader if phase == 'train' else valid_dataloader

#             for inputs, labels in dataloader:
#                 inputs = inputs.to(device)
#                 labels = labels.to(device).float().unsqueeze(1)

#                 optimizer.zero_grad()

#                 with torch.set_grad_enabled(phase == 'train'):
#                     outputs = model(inputs)
#                     preds = torch.sigmoid(outputs) >= 0.5
#                     loss = criterion(outputs, labels)

#                     if phase == 'train':
#                         loss.backward()
#                         optimizer.step()

#                 running_loss += loss.item() * inputs.size(0)
#                 running_corrects += torch.sum(preds == labels.data)

#             epoch_loss = running_loss / len(dataloader.dataset)
#             epoch_acc = running_corrects.double() / len(dataloader.dataset)

#             print(f'{phase} Loss: {epoch_loss:.4f} Acc: {epoch_acc:.4f}')

#             if phase == 'val' and epoch_acc > best_acc:
#                 best_acc = epoch_acc
#                 best_model_wts = copy.deepcopy(model.state_dict())

#     print(f'Best val Acc: {best_acc:.4f}')
#     model.load_state_dict(best_model_wts)
#     return model

# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# model = model.to(device)


# # Training
# num_epochs = 10
# for epoch in range(num_epochs):
#     model.train()
#     epoch_loss = 0
#     for images, labels in train_dataloader:
#         images, labels = images.to(device), labels.to(device).float().unsqueeze(1)

#         optimizer.zero_grad()
#         outputs = model(images)
#         loss = criterion(outputs, labels)
#         loss.backward()
#         optimizer.step()    

#%%

# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# model = model.to(device)


import torch.nn.functional as F

def train_model(model, train_dataloader, valid_dataloader, criterion, optimizer, num_epochs=25):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    for epoch in range(num_epochs):
        model.train()  # Set model to training mode

        running_loss = 0.0
        running_corrects = 0
        all_labels = []
        all_preds = []

        # Iterate over data.
        for inputs, labels in train_dataloader:
            inputs = inputs.to(device)
            labels = labels.to(device).float()

            # Zero the parameter gradients
            optimizer.zero_grad()

            # Forward
            outputs = model(inputs)
            # outputs = outputs.squeeze()
            if outputs.shape == torch.Size([1, 1]):  # outputs is a scalar           
                outputs = outputs.squeeze().unsqueeze(0)   # Reshape scalar to [1]                
            else:
                outputs = outputs.squeeze()    # Squeeze the output if necessary

            loss = criterion(outputs, labels)
            preds = outputs.sigmoid() > 0.5

            # Backward + optimize
            loss.backward()
            optimizer.step()

            # Statistics
            running_loss += loss.item() * inputs.size(0)
            running_corrects += torch.sum(preds == labels.data)
            all_labels.append(labels.cpu().numpy())
            all_preds.append(preds.cpu().numpy())

        epoch_loss = running_loss / len(train_dataloader.dataset)
        epoch_acc = running_corrects.double() / len(train_dataloader.dataset)

        all_labels = np.concatenate(all_labels)
        all_preds = np.concatenate(all_preds)
        cm = confusion_matrix(all_labels, all_preds)
        tn, fp, fn, tp = cm.ravel()
        class_0_accuracy = tn / (tn + fp)
        class_1_accuracy = tp / (tp + fn)

        print(f'Epoch {epoch + 1}/{num_epochs} - Training Loss: {epoch_loss:.4f} Acc: {epoch_acc:.4f}')
        print(f'Class 0 Accuracy: {class_0_accuracy:.4f}, Class 1 Accuracy: {class_1_accuracy:.4f}')

        # Validation phase
        model.eval()  # Set model to evaluate mode
        valid_loss = 0.0
        valid_corrects = 0
        valid_labels = []
        valid_preds = []

        for inputs, labels in valid_dataloader:
            inputs = inputs.to(device)
            labels = labels.to(device).float()

            with torch.no_grad():
                outputs = model(inputs)
                # outputs = outputs.squeeze()
                
                if outputs.shape == torch.Size([1, 1]):  # outputs is a scalar
                    outputs = outputs.squeeze().unsqueeze(0)   # Reshape scalar to [1]                
                else:
                    outputs = outputs.squeeze()    # Squeeze the output if necessary

                loss = criterion(outputs, labels)
                preds = outputs.sigmoid() > 0.5

            valid_loss += loss.item() * inputs.size(0)
            valid_corrects += torch.sum(preds == labels.data)
            valid_labels.append(labels.cpu().numpy())
            valid_preds.append(preds.cpu().numpy())

        valid_epoch_loss = valid_loss / len(valid_dataloader.dataset)
        valid_epoch_acc = valid_corrects.double() / len(valid_dataloader.dataset)

        valid_labels = np.concatenate(valid_labels)
        valid_preds = np.concatenate(valid_preds)
        cm = confusion_matrix(valid_labels, valid_preds)
        tn, fp, fn, tp = cm.ravel()
        valid_class_0_accuracy = tn / (tn + fp)
        valid_class_1_accuracy = tp / (tp + fn)

        print(f'Epoch {epoch + 1}/{num_epochs} - Validation Loss: {valid_epoch_loss:.4f} Acc: {valid_epoch_acc:.4f}')
        print(f'Validation Class 0 Accuracy: {valid_class_0_accuracy:.4f}, Validation Class 1 Accuracy: {valid_class_1_accuracy:.4f}')

    return model





#%%


# Call the training function
trained_model = train_model(model, train_dataloader, valid_dataloader, criterion, optimizer, num_epochs=num_epochs)


#%%

# for inputs, labels in train_dataloader:
#     print(inputs)
#     print(labels)
