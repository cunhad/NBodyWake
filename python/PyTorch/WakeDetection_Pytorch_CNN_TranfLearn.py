#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:41:50 2024

@author: Disrael
"""

#%%

# paremeters

batch_size = 2
num_workers = 3
TRAIN_RATIO = 0.8       #fraction of total dataset that will go to train+validation
VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation
OUTPUT_DIM = 1          # 2 classes for classification labels
SEED = 1234
pretrained_size = 512


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
# import sys
# path_analy = os.getcwd() +'/' 
# sys.path.append(path_analy)
from DataCleaning import file_list,list_all_files,balanced_list_of_files

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
            transform (callable, optional): Optional transform to be applied
                on a sample.
        """
        self.file_list = file_list
        self.transform = transform

        # Extract classes from folder structure
        self.classes = list(set([os.path.basename(os.path.dirname(file_path)) for file_path in file_list]))
        self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}
        self.labels = [self.class_to_idx[os.path.basename(os.path.dirname(file_path))] for file_path in file_list]


    def __len__(self):
        return len(self.file_list)

    def __getitem__(self, idx):
        img_path = self.file_list[idx]
        image = Image.open(img_path).convert("RGB")
        class_name = os.path.basename(os.path.dirname(img_path))
        label = self.class_to_idx[class_name]

        if self.transform:
            image = self.transform(image)

        return image, label

# Load datasets

folder_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50_test/"
files_list = list_all_files(folder_path)

valid_data_, test_train_data_,selected_indices_val,selected_indices_test_train = file_list(folder_path,VALID_RATIO)


folder_path_balance = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr50/"
files_list_balance  = list_all_files(folder_path_balance)


files_list_balanced_val = balanced_list_of_files(files_list,files_list_balance,valid_data_,selected_indices_val)
files_list_balanced_test_train = balanced_list_of_files(files_list,files_list_balance,test_train_data_,selected_indices_test_train)


# valid_data_, test_train_data_,_ = file_list(folder_path,VALID_RATIO)
valid_data__ = CustomImageDataset(file_list=files_list_balanced_val, transform=test_transforms)
test_train_data__ = CustomImageDataset(file_list=files_list_balanced_test_train, transform=test_transforms)


# valid_data_, test_train_data_,_ = file_list(folder_path,VALID_RATIO)
# valid_data__ = CustomImageDataset(file_list=valid_data_, transform=test_transforms)
# test_train_data__ = CustomImageDataset(file_list=test_train_data_, transform=test_transforms)


n_train_examples = int(len(test_train_data__) * TRAIN_RATIO)
n_test_examples = len(test_train_data__) - n_train_examples

train_data_, test_data_ = data.random_split(test_train_data__,
                                           [n_train_examples, n_test_examples])

# Create data loaders.

valid_dataloader = DataLoader(valid_data__, batch_size=batch_size, shuffle=False, num_workers=num_workers)
test_dataloader = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
train_dataloader = DataLoader(train_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)


# # Create data loaders.
# train_dataloader = DataLoader(train_data, batch_size=batch_size, shuffle=True, num_workers=num_workers)
# test_dataloader = DataLoader(test_data, batch_size=batch_size, shuffle=False, num_workers=num_workers)
# valid_dataloader = DataLoader(valid_data, batch_size=batch_size, shuffle=False, num_workers=num_workers)

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




#%%

print(test_train_data__.class_to_idx)   # this is the label dictionary



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

# Print the model architecture
print(model)

# Freeze all layers
for param in model.parameters():
    param.requires_grad = False
    
    
# Modify the final layer for binary classification
IN_FEATURES = model.classifier[-1].in_features
model.classifier[-1] = nn.Linear(IN_FEATURES, OUTPUT_DIM)   

# Unfreeze the last layer
for param in model.classifier[-1].parameters():
    param.requires_grad = True 
    
# Define the loss function and optimizer
criterion = nn.BCEWithLogitsLoss()  # Binary Cross Entropy with Logits Loss
optimizer = optim.Adam(model.classifier[-1].parameters(), lr=0.001)


#%%

# Print the model architecture
print(model)    


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

# def train_model(model, train_dataloader, valid_dataloader, criterion, optimizer, num_epochs=25):
#     device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#     model.to(device)

#     for epoch in range(num_epochs):
#         model.train()  # Set model to training mode

#         running_loss = 0.0
#         running_corrects = 0

#         # Iterate over data.
#         for inputs, labels in train_dataloader:
#             inputs = inputs.to(device)
#             labels = labels.to(device).float()
#             # labels = labels.to(device).float().view(-1, 1)  # Ensure labels are the correct shape
#             # labels = labels.to(device).float().view(-1, 1).squeeze()  # Ensure labels are the correct shape

                     

#             # Zero the parameter gradients
#             optimizer.zero_grad()

#             # Forward
#             outputs = model(inputs)
            
#             # print(outputs.shape)
            
#             if outputs.shape == torch.Size([1, 1]):  # outputs is a scalar
#                 # print("this")
#                 # outputs = outputs.squeeze()    # Reshape scalar to [1]                
#                 outputs = outputs.squeeze().unsqueeze(0)   # Reshape scalar to [1]                
#             else:
#                 # print("that")
#                 # outputs = torch.tensor([outputs.item()], device='cuda:0')
#                 outputs = outputs.squeeze()    # Squeeze the output if necessary
            
#             # if outputs.dim() != 0:  # outputs is a scalar
#             #     outputs = outputs.squeeze()  # Squeeze the output if necessary
            
#             # outputs = outputs.squeeze()  # Squeeze the output if necessary
#             # print(labels)   
#             # print(outputs) 
            
#             # outputs = outputs.squeeze()  # Squeeze the output if necessary
#             # outputs = outputs.view(-1, 1).squeeze()   # Ensure outputs are also [batch_size, 1] if not already
#             loss = criterion(outputs, labels)
#             preds = outputs.sigmoid() > 0.5
#             # Backward + optimize
#             loss.backward()
#             optimizer.step()

#             # Statistics
#             running_loss += loss.item() * inputs.size(0)
#             running_corrects += torch.sum(preds == labels.data)

#         epoch_loss = running_loss / len(train_dataloader.dataset)
#         epoch_acc = running_corrects.double() / len(train_dataloader.dataset)

#         print(f'Epoch {epoch + 1}/{num_epochs} - Training Loss: {epoch_loss:.4f} Acc: {epoch_acc:.4f}')

#         # Validation phase
#         model.eval()  # Set model to evaluate mode
#         valid_loss = 0.0
#         valid_corrects = 0

#         for inputs, labels in valid_dataloader:
#             inputs = inputs.to(device)
#             # labels = labels.to(device)
#             labels = labels.to(device).float()

#             with torch.no_grad():
#                 outputs = model(inputs)
                
#                 if outputs.shape == torch.Size([1, 1]):  # outputs is a scalar
#                     # print("this")
#                     # outputs = outputs.squeeze()    # Reshape scalar to [1]                
#                     outputs = outputs.squeeze().unsqueeze(0)   # Reshape scalar to [1]                
#                 else:
#                     # print("that")
#                     # outputs = torch.tensor([outputs.item()], device='cuda:0')
#                     outputs = outputs.squeeze()    # Squeeze the output if necessary
                
#                 loss = criterion(outputs, labels)

#                 preds = outputs.sigmoid() > 0.5

#             valid_loss += loss.item() * inputs.size(0)
#             valid_corrects += torch.sum(preds == labels.data)

#         valid_epoch_loss = valid_loss / len(valid_dataloader.dataset)
#         valid_epoch_acc = valid_corrects.double() / len(valid_dataloader.dataset)

#         print(f'Epoch {epoch + 1}/{num_epochs} - Validation Loss: {valid_epoch_loss:.4f} Acc: {valid_epoch_acc:.4f}')

#     return model

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
trained_model = train_model(model, train_dataloader, valid_dataloader, criterion, optimizer, num_epochs=10)


#%%

# for inputs, labels in train_dataloader:
#     print(inputs)
#     print(labels)
