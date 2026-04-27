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
import ast


parser = argparse.ArgumentParser(description='efficientnet_b7 wake classification model')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--path_data', type=str, help='')
parser.add_argument('--path_void', type=str, help='')
parser.add_argument('--path_WakeSignal', type=str, help='')

parser.add_argument('--n_angle', type=int, default=96, help='')
parser.add_argument('--rangeSampl', type=parse_range, default='5001-5100', help="Range of samples ids (e.g., '5001-5100,5200,5211')")

parser.add_argument('--slices_void', type=int, default=33, help='')
parser.add_argument('--slices_signal', type=int, default=8, help='')
# parser.add_argument('--percentage_positiveWakeSig', type=float, default=10, help='')
parser.add_argument('--validation_fraction', type=float, default=0.1, help='')
parser.add_argument('--train_tt_fraction', type=float, default=0.8, help='')
parser.add_argument('--wake_top_percentage', type=float, default=10, help='')
parser.add_argument('--void_percentage', type=float, default=0, help='')


parser.add_argument('--batch_size', type=int, default=4, help='')
parser.add_argument('--num_workers', type=int, default=0, help='')
parser.add_argument('--num_epochs', type=int, default=50, help='')

parser.add_argument('--Nmesh', default='[512,512,32]', help='')


args = parser.parse_args()


# parameters

path_data = args.path_data
print("File Path in = "+ str(path_data))

path_void = args.path_void
print("File Path Void in = "+ str(path_void))

path_WakeSignal = args.path_WakeSignal
print("File Path Wake Signal in = "+ str(path_WakeSignal))


n_angle = args.n_angle
print("angles each sample = "+ str(n_angle))

# Access the parsed values
rang = args.rangeSampl
print("Samples = ", rang)




slices_void = args.slices_void
print("Slices for void = "+ str(slices_void))

slices_signal = args.slices_signal
print("Slices for wake signal = ", slices_signal)

validation_fraction = args.validation_fraction
print("Validation fraction of total data = ", validation_fraction)

train_tt_fraction = args.train_tt_fraction
print("Train fraction of train + test data = ", train_tt_fraction)


wake_top_percentage = args.wake_top_percentage
print("Percentage top signal wake = ", wake_top_percentage)

void_percentage = args.void_percentage
print("Percentage lower voids wake = ", void_percentage)




batch_size =  args.batch_size
print("Batch size = "+ str(batch_size))

num_workers = args.num_workers
print("Num of CPU workers = "+ str(num_workers))

num_epochs = args.num_epochs
print("Num epochs = "+ str(num_epochs))

Nmesh =  ast.literal_eval(args.Nmesh)
print("Nmesh= "+ str(Nmesh))




# # parameters

# path_data = args.path_data
# path_data = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/"
# print("File Path in = "+ str(path_data))

# path_void = args.path_void
# path_void = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat/void/"
# print("File Path Void in = "+ str(path_void))

# path_WakeSignal = args.path_WakeSignal
# path_WakeSignal = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_stat_2dc1l1_3dc1l1/"
# print("File Path Wake Signal in = "+ str(path_WakeSignal))




# n_angle = args.n_angle
# # num_epochs = 1
# print("angles each sample = "+ str(n_angle))

# # Access the parsed values
# rang = args.rangeSampl
# # rangeAng = parse_range('5001-5101')
# print("Samples = ", rang)




# slices_void = args.slices_void
# # slices_signal = 33
# print("Slices for void = "+ str(slices_void))

# slices_signal = args.slices_signal
# # slices_signal = 32
# print("Slices for wake signal = ", slices_signal)

# # percentage_positiveWakeSig = args.percentage_positiveWakeSig
# # # percentage_positiveWakeSig = 10
# # print("Percentage positive Wake Signal = ", percentage_positiveWakeSig)

# validation_fraction = args.validation_fraction
# # validation_fraction = 0.1
# print("Validation fraction of total data = ", validation_fraction)

# train_tt_fraction = args.train_tt_fraction
# # train_tt_fraction = 0.1
# print("Train fraction of train + test data = ", train_tt_fraction)


# wake_top_percentage = args.wake_top_percentage
# # wake_top_percentage = 10
# print("Percentage top signal wake = ", wake_top_percentage)

# void_percentage = args.void_percentage
# # void_percentage = 0.1
# print("Percentage lower voids wake = ", void_percentage)




# batch_size =  args.batch_size
# batch_size =  1
# print("Batch size = "+ str(batch_size))

# num_workers = args.num_workers
# # num_workers = 0
# print("Num of CPU workers = "+ str(num_workers))

# num_epochs = args.num_epochs
# # num_epochs = 1
# print("Num epochs = "+ str(num_epochs))

# Nmesh =  ast.literal_eval(args.Nmesh)
# # Nmesh = [48,48,48]
# # Nmesh = [48,48,12]
# # Nmesh = [512,512,32]
# # Nmesh = [512,512,512]
# print("Nmesh= "+ str(Nmesh))












# # batch_size = 32
# TRAIN_RATIO = 0.8       #fraction of total dataset that will go to train+validation
# VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation
OUTPUT_DIM = 1          # 2 classes for classification labels
# SEED = 1234
pretrained_size = 512
# pretrained_size = 160


wake_spec = ["4Mpc_2048c_1024p_zi63_nowakem/","4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"]

WAKE_NAME   = "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m"
NOWAKE_NAME = "4Mpc_2048c_1024p_zi63_nowakem"









#%%

# Import the necessary libraries:

# Torch stuff
    
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, random_split, Subset
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

from scipy.ndimage import gaussian_filter


# Import the spliter function from DataCleaning.py'
import sys

# from DataCleaning import file_list,list_all_files,balanced_list_of_files
path_analy = os.getcwd() +'/' 
# path_analy = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path_analy)
from Data_preSelection  import  *

# keep track of true positives, false positives, true negatives, and false negatives for each class
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit


# To count the number of images in each label
from collections import Counter




#%%

# prioritize wake recall
POS_WEIGHT = 0.7  # try 0.7, 0.8, 1.0 — never below 0.5
# THRESH_GRID = np.linspace(0.05, 0.95, 37)
THRESH_GRID = np.arange(0.30, 0.80, 0.01)
DEFAULT_THRESHOLD = 0.50
FULL_VAL_EVERY = 5   # evaluate full samples_val every 5 epochs

#%%

# Define the data transformations:

    
# pretrained_means = [0.485, 0.456, 0.406]
# pretrained_stds = [0.229, 0.224, 0.225]

# train_transforms = transforms.Compose([
#                            transforms.Resize(pretrained_size),
#                            # transforms.RandomRotation(5),
#                            transforms.RandomHorizontalFlip(0.5),
#                            # transforms.RandomCrop(pretrained_size, padding=10),
#                            transforms.ToTensor(),
#                            transforms.Normalize(mean=pretrained_means,
#                                                 std=pretrained_stds)
#                        ])

# test_transforms = transforms.Compose([
#                            transforms.Resize(pretrained_size),
#                            transforms.ToTensor(),
#                            transforms.Normalize(mean=pretrained_means,
#                                                 std=pretrained_stds)
#                        ])



#%%

def load_bin_volume(
    filename: str,
    grid_shape: tuple[int, int, int],
    dtype: np.dtype = np.float32,
    apply_matlab_transform: bool = True,
) -> np.ndarray:
    """
    Lê o volume 3D completo de um .bin, reshape em (D,H,W),
    e opcionalmente aplica a transformação estilo Matlab:

        dc  = (vol - mean_global) / mean_global
        vol = atan((dc + 1) * 16)

    Retorna: vol com shape (D, H, W), dtype float32.
    """
    Nx, Ny, Nz = grid_shape
    n_elems = Nx * Ny * Nz

    flat = np.fromfile(filename, dtype=dtype, count=n_elems)
    if flat.size != n_elems:
        raise ValueError(
            f"Arquivo {filename} tem {flat.size} elementos, "
            f"mas esperados {n_elems} para grid_shape={grid_shape}."
        )

    # (Nz, Nx, Ny) assumindo que cada slice era Nx x Ny em C-order
    vol = flat.reshape(Nz, Nx, Ny)
    # Transpõe para (D=Nz, H=Ny, W=Nx), consistente com o .T que você usava para slices
    vol = vol.transpose(0, 2, 1)  # (D, H, W)

    vol = vol.astype(np.float32)

    if apply_matlab_transform:
        mean_global = float(vol.mean())
        if mean_global != 0.0:
            dc  = (vol - mean_global) / mean_global
            vol = np.arctan((dc + 1.0) * 16.0).astype(np.float32)
            
    # #only for debug
    # vol = vol[:,::4,::4]

    return vol


## verify if it loads correctly


class BinVolumeDataset(Dataset):
    """
    Cada item = volume 3D completo de um .bin.
    Saída: (tensor (1, D, H, W), label 0/1).
    """

    def __init__(
        self,
        bin_paths: list[str],
        wake_infos: list[str],
        grid_shape: tuple[int, int, int],
        subvol_labels=None,
        dtype: np.dtype = np.float32,
    ):
        assert len(bin_paths) == len(wake_infos), \
            "bin_paths e wake_infos devem ter o mesmo comprimento."

        self.bin_paths  = bin_paths
        self.wake_infos = wake_infos
        self.grid_shape = tuple(grid_shape)
        self.dtype      = dtype
        self.subvol_labels = subvol_labels

        if self.subvol_labels is not None:
            assert len(self.subvol_labels) == len(self.bin_paths), \
                "subvol_labels must have same length as bin_paths."

        self.class_to_idx = {
            NOWAKE_NAME: 0,
            WAKE_NAME:   1,
        }

        self.classes = [None] * len(self.class_to_idx)
        for name, idx in self.class_to_idx.items():
            self.classes[idx] = name

        self.labels = [self._label_from_wake_info(wi) for wi in wake_infos]

    def _label_from_wake_info(self, wi: str) -> int:
        if WAKE_NAME in wi:
            return 1
        elif NOWAKE_NAME in wi:
            return 0
        else:
            raise ValueError(f"wake_info desconhecido: {wi}")

    def __len__(self):
        return len(self.bin_paths)

    def __getitem__(self, idx):
        bin_path = self.bin_paths[idx]
        label    = self.labels[idx]

        # 1) lê volume e aplica a transformação Matlab internamente
        vol = load_bin_volume(
            filename           = bin_path,
            grid_shape         = self.grid_shape,
            dtype              = self.dtype,
            apply_matlab_transform = True,
        )  # (D, H, W)

        # vol = (vol - vol.mean()) / (vol.std() + 1e-6)
        # # # 2) normalização opcional por volume (exemplo: min–max)
        # # vmin, vmax = vol.min(), vol.max()
        # # if vmax > vmin:
        # #     vol_norm = (vol - vmin) / (vmax - vmin)
        # # else:
        # #     vol_norm = np.zeros_like(vol, dtype=np.float32)

        # # 3) tensor (1, D, H, W)
        # # vol_tensor = torch.from_numpy(vol_norm).unsqueeze(0)
        # vol_tensor = torch.from_numpy(vol).unsqueeze(0)

        # channel 1: robust z-score
        vol = (vol - vol.mean()) / (vol.std() + 1e-6)
        vol = np.clip(vol, -5.0, 5.0)
        
        # channel 2: high-pass residual
        smooth = gaussian_filter(vol, sigma=(0.7, 1.2, 1.2))
        high = vol - smooth
        high = np.clip(high, -5.0, 5.0)
        
        # vol_tensor = torch.from_numpy(np.stack([vol, high], axis=0)).float()
        
        # final tensor: (2, D, H, W)
        vol_2ch = np.stack([vol, high], axis=0).astype(np.float32, copy=False)
        vol_tensor = torch.from_numpy(np.ascontiguousarray(vol_2ch))


        # return vol_tensor, label

        vol_label = torch.tensor(label, dtype=torch.float32)

        if self.subvol_labels is None:
            subvol_label = torch.zeros(slices_signal, dtype=torch.float32)
        else:
            subvol_label = torch.tensor(self.subvol_labels[idx], dtype=torch.float32)
        
        return vol_tensor, vol_label, subvol_label



#%%


# Load datasets

files_list_all = list_all_files(path_data, rang)
samples_all, anglids_all, wake_infos_all = extract_info_bin(files_list_all)

# only keep files in rang

rang_set = set(rang)

filtered = [
    (f, s, a, w)
    for f, s, a, w in zip(files_list_all, samples_all, anglids_all, wake_infos_all)
    if s in rang_set
]

files_list_all = [x[0] for x in filtered]
samples_all    = [x[1] for x in filtered]
anglids_all    = [x[2] for x in filtered]
wake_infos_all = [x[3] for x in filtered]



all_data_nowake_void, all_data_wake_void = data_void_out3d(rang, n_angle, slices_void, wake_spec, path_void)
all_data_nowake_signal, all_data_wake_signal = data_signal_out3d(rang, n_angle, slices_signal, wake_spec, path_WakeSignal)

all_data_nowake_signal_subvol, all_data_wake_signal_subvol = data_signal_subvol_out3d(
    rang, n_angle, slices_signal, wake_spec, path_WakeSignal
)



data_void = extract_stat3d(samples_all, anglids_all, wake_infos_all, all_data_nowake_void, all_data_wake_void, rang)
data_signal = extract_stat3d(samples_all, anglids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, rang)
data_signal_diff = extract_stat_diff3d(samples_all, anglids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, rang)


#%%


# split val and trainTest

# first split per sample (so the validation is independent)

validation_indices, train_test_indices = split_unique_samples(samples_all, wake_infos_all, validation_fraction)

files_list_validation = [files_list_all[i] for i in validation_indices]
files_list_trainTest = [files_list_all[i] for i in train_test_indices]

samples_val    = [samples_all[i] for i in validation_indices]
anglids_val    = [anglids_all[i] for i in validation_indices]
wake_infos_val = [wake_infos_all[i] for i in validation_indices]
data_void_val =  [data_void[i] for i in validation_indices]
data_signal_val =  [data_signal[i] for i in validation_indices]
data_signal_diff_val = [data_signal_diff[i] for i in validation_indices]

samples_tt    = [samples_all[i] for i in train_test_indices]
anglids_tt    = [anglids_all[i] for i in train_test_indices]
wake_infos_tt = [wake_infos_all[i] for i in train_test_indices]
data_void_tt    = [data_void[i] for i in train_test_indices]
data_signal_tt = [data_signal[i] for i in train_test_indices]
data_signal_diff_tt = [data_signal_diff[i] for i in train_test_indices]

#%%




# validation

# Select top 50% wake values, keep void_percentage 
selected_files_val, selected_positions_val, selected_signal_diff_val = select_extreme_files(
    data_signal_diff_val, wake_infos_val, files_list_validation, data_void_val,
    wake_top_percentage, void_percentage
)

find_extreme_files(selected_signal_diff_val, selected_files_val)
selected_wake_infos_val = [wake_infos_val[i] for i in selected_positions_val]



#%%

# for subvolumes

selected_samples_val = [samples_val[i] for i in selected_positions_val]
selected_anglids_val = [anglids_val[i] for i in selected_positions_val]

def build_subvol_labels(
    samples,
    anglids,
    wake_infos,
    all_data_nowake_signal_subvol,
    all_data_wake_signal_subvol,
    rang,
    slices_signal,
    eps=1e-6,
):
    sample_to_idx = {s: i for i, s in enumerate(rang)}
    subvol_labels = []

    for sample, anglid, wake_info in zip(samples, anglids, wake_infos):
        sample_id = sample_to_idx[sample]
        angle_id = anglid - 1

        if wake_info == WAKE_NAME:
            diff = (
                all_data_wake_signal_subvol[sample_id, angle_id, :]
                - all_data_nowake_signal_subvol[sample_id, angle_id, :]
            )

            diff = np.nan_to_num(diff, nan=0.0)
            diff = np.maximum(diff, 0.0)

            maxval = diff.max()
            if maxval > eps:
                label8 = diff / maxval
            else:
                label8 = np.zeros(slices_signal, dtype=np.float32)

        elif wake_info == NOWAKE_NAME:
            label8 = np.zeros(slices_signal, dtype=np.float32)

        else:
            raise ValueError(f"Unknown wake_info: {wake_info}")

        subvol_labels.append(label8.astype(np.float32))

    return subvol_labels

selected_subvol_labels_val = build_subvol_labels(
    selected_samples_val,
    selected_anglids_val,
    selected_wake_infos_val,
    all_data_nowake_signal_subvol,
    all_data_wake_signal_subvol,
    rang,
    slices_signal,
)





#%%




# train and test

# Select top 50% wake values, keep void_percentage 
selected_files_tt, selected_positions_tt, selected_signal_diff_tt = select_extreme_files(
    data_signal_diff_tt, wake_infos_tt, files_list_trainTest, data_void_tt,
    wake_top_percentage, void_percentage
)

# selected_sliceids_tt   = [sliceids_tt[i] for i in selected_positions_tt]
selected_wake_infos_tt = [wake_infos_tt[i] for i in selected_positions_tt]


find_extreme_files(selected_signal_diff_tt, selected_files_tt)

#%%

# for subvoumes:

selected_samples_tt = [samples_tt[i] for i in selected_positions_tt]
selected_anglids_tt = [anglids_tt[i] for i in selected_positions_tt]

selected_subvol_labels_tt = build_subvol_labels(
    selected_samples_tt,
    selected_anglids_tt,
    selected_wake_infos_tt,
    all_data_nowake_signal_subvol,
    all_data_wake_signal_subvol,
    rang,
    slices_signal,
)


#%%

# shape of your 3D grid in each bin file:
# grid_shape = (slices_signal, original_H, original_W)
# e.g. if you know it's 32x512x512:
# grid_shape = (slices_signal, 512, 512)

# valid_data__ = BinVolumeDataset(
#     bin_paths=selected_files_val,
#     wake_infos=selected_wake_infos_val,
#     grid_shape=Nmesh,
# )

# test_train_data__ = BinVolumeDataset(
#     bin_paths=selected_files_tt,
#     wake_infos=selected_wake_infos_tt,
#     grid_shape=Nmesh,
# )

valid_data__ = BinVolumeDataset(
    bin_paths=selected_files_val,
    wake_infos=selected_wake_infos_val,
    grid_shape=Nmesh,
    subvol_labels=selected_subvol_labels_val,
)

test_train_data__ = BinVolumeDataset(
    bin_paths=selected_files_tt,
    wake_infos=selected_wake_infos_tt,
    grid_shape=Nmesh,
    subvol_labels=selected_subvol_labels_tt,
)

#%%



# this should be uncommented (if not on debug)


n_train_examples = min(int(len(test_train_data__) * train_tt_fraction),int(len(test_train_data__)) - 2)
n_test_examples = len(test_train_data__) - n_train_examples

labels_tt = test_train_data__.labels
sss = StratifiedShuffleSplit(
    n_splits=1,
    test_size=n_test_examples,
    train_size=n_train_examples,
    random_state=42
)

train_idx, test_idx = next(sss.split(range(len(labels_tt)), labels_tt))

train_data_ = Subset(test_train_data__, train_idx)
test_data_  = Subset(test_train_data__, test_idx)



# # this should be uncommented (if on debug)


# print("\n=== TINY OVERFIT TEST ===")

# labels_tt = test_train_data__.labels

# idx_wake = [i for i, y in enumerate(labels_tt) if y == 1][:4]
# idx_nowake = [i for i, y in enumerate(labels_tt) if y == 0][:4]

# tiny_indices = idx_wake + idx_nowake

# print("Tiny wake indices:", idx_wake)
# print("Tiny no-wake indices:", idx_nowake)
# print("Tiny dataset size:", len(tiny_indices))

# train_data_ = Subset(test_train_data__, tiny_indices)
# test_data_  = Subset(test_train_data__, tiny_indices)


# # # medium debug

# # after you create train_data_
# n_probe = 64
# probe_idx = np.random.RandomState(42).choice(len(train_data_), size=min(n_probe, len(train_data_)), replace=False)
# probe_data = Subset(train_data_, probe_idx)
# probe_loader = DataLoader(probe_data, batch_size=1, shuffle=True, num_workers=num_workers)

# # temporarily train only on probe_loader
# train_dataloader = probe_loader
# valid_dataloader = probe_loader



#%%


# Create data loaders.

# this should be uncommented (if not on debug)


valid_dataloader = DataLoader(valid_data__, batch_size=batch_size, shuffle=False, num_workers=num_workers)
test_dataloader = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
train_dataloader = DataLoader(train_data_, batch_size=batch_size, shuffle=True, num_workers=num_workers)


# # this should be uncommented (if on debug)


# valid_dataloader = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
# test_dataloader  = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
# train_dataloader = DataLoader(train_data_, batch_size=batch_size, shuffle=True,  num_workers=num_workers)

for X, vol_label, subvol_label in train_dataloader:
    print(f"Shape of X [B, C, D, H, W]: {X.shape}")
    print(f"Shape of vol_label: {vol_label.shape} {vol_label.dtype}")
    print(f"Shape of subvol_label: {subvol_label.shape} {subvol_label.dtype}")
    break
    
for X, vol_label, subvol_label  in test_dataloader:
    print(f"Shape of X [B, C, D, H, W]: {X.shape}")
    print(f"Shape of vol_label: {vol_label.shape} {vol_label.dtype}")
    print(f"Shape of subvol_label: {subvol_label.shape} {subvol_label.dtype}")
    break  

for X, vol_label, subvol_label  in valid_dataloader:
    print(f"Shape of X [B, C, D, H, W]: {X.shape}")
    print(f"Shape of vol_label: {vol_label.shape} {vol_label.dtype}")
    print(f"Shape of subvol_label: {subvol_label.shape} {subvol_label.dtype}")
    break  



print(test_train_data__.class_to_idx)   # this is the label dictionary









#%%


def count_class_volumes(subset, original_dataset):
    labels = [original_dataset.labels[i] for i in subset.indices]
    class_counts = Counter(labels)
    class_names = original_dataset.classes
    for class_idx, count in class_counts.items():
        class_name = class_names[class_idx]
        print(f"Class: {class_name}, Number of volumes: {count}")

print("\nTraining Data:")
count_class_volumes(train_data_, test_train_data__)

print("\nTesting Data:")
count_class_volumes(test_data_, test_train_data__)



# Count the number of images in each class
class_counts = Counter(valid_data__.labels)
# Get the class names
class_names = valid_data__.classes

# Display the number of images in each class
print("\nValidation Data:")
for class_idx, count in class_counts.items():
    class_name = class_names[class_idx]
    print(f"Class: {class_name}, Number of volumes: {count}")


#%%


# class DoubleConv3D(nn.Module):
#     def __init__(self, in_channels, out_channels, norm="group", num_groups=4):
#         super().__init__()

#         def make_norm(num_channels):
#             if norm == "group":
#                 # num_groups must divide num_channels
#                 groups = min(num_groups, num_channels)
#                 while num_channels % groups != 0 and groups > 1:
#                     groups -= 1
#                 return nn.GroupNorm(num_groups=groups, num_channels=num_channels)

#             elif norm == "instance":
#                 return nn.InstanceNorm3d(num_channels)

#             elif norm == "batch":
#                 return nn.BatchNorm3d(num_channels)

#             else:
#                 raise ValueError("norm must be 'group', 'instance', or 'batch'")

#         self.block = nn.Sequential(
#             nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1),
#             make_norm(out_channels),
#             nn.ReLU(inplace=True),

#             nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1),
#             make_norm(out_channels),
#             nn.ReLU(inplace=True),
#         )

#     def forward(self, x):
#         return self.block(x)


# class UNet3DClassifier(nn.Module):
#     """
#     3D U-Net-style encoder classifier.
#     Input:  (B, 1, D, H, W)
#     Output: (B,) logits
#     """

#     def __init__(self, in_channels=1, base_channels=16, norm="instance"):
#         super().__init__()

#         # Encoder
#         self.enc1 = DoubleConv3D(in_channels, base_channels, norm=norm)
#         self.pool1 = nn.MaxPool3d(kernel_size=(2, 2, 2))

#         self.enc2 = DoubleConv3D(base_channels, base_channels * 2, norm=norm)
#         self.pool2 = nn.MaxPool3d(kernel_size=(2, 2, 2))

#         self.enc3 = DoubleConv3D(base_channels * 2, base_channels * 4, norm=norm)
#         self.pool3 = nn.MaxPool3d(kernel_size=(2, 2, 2))

#         # Bottleneck
#         self.bottleneck = DoubleConv3D(base_channels * 4, base_channels * 8, norm=norm)

#         # Classification head
#         self.global_pool = nn.AdaptiveAvgPool3d(1)   # -> (B, C, 1,1,1)
#         self.classifier = nn.Linear(base_channels * 8, 1)

#     def forward(self, x):
#         x = self.enc1(x)      # (B,16,32,512,512)
#         x = self.pool1(x)     # (B,16,16,256,256)

#         x = self.enc2(x)      # (B,32,16,256,256)
#         x = self.pool2(x)     # (B,32, 8,128,128)

#         x = self.enc3(x)      # (B,64, 8,128,128)
#         x = self.pool3(x)     # (B,64, 4, 64, 64)

#         x = self.bottleneck(x)  # (B,128,4,64,64)

#         x = self.global_pool(x)   # (B,128,1,1,1)
#         x = x.view(x.size(0), -1) # (B,128)
#         x = self.classifier(x)    # (B,1)
#         x = x.squeeze(1)          # (B,)

#         return x

#%%
    
# class ResidualBlock3D(nn.Module):
#     def __init__(self, in_channels, out_channels, norm="group", num_groups=4):
#         super().__init__()

#         def make_norm(num_channels):
#             if norm == "group":
#                 groups = min(num_groups, num_channels)
#                 while num_channels % groups != 0 and groups > 1:
#                     groups -= 1
#                 return nn.GroupNorm(num_groups=groups, num_channels=num_channels)
#             elif norm == "instance":
#                 return nn.InstanceNorm3d(num_channels)
#             elif norm == "batch":
#                 return nn.BatchNorm3d(num_channels)
#             else:
#                 raise ValueError("norm must be 'group', 'instance', or 'batch'")

#         self.conv1 = nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1, bias=False)
#         self.norm1 = make_norm(out_channels)
#         self.relu1 = nn.ReLU(inplace=True)

#         self.conv2 = nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1, bias=False)
#         self.norm2 = make_norm(out_channels)

#         if in_channels != out_channels:
#             self.skip = nn.Conv3d(in_channels, out_channels, kernel_size=1, bias=False)
#         else:
#             self.skip = nn.Identity()

#         self.relu2 = nn.ReLU(inplace=True)

#     def forward(self, x):
#         identity = self.skip(x)

#         out = self.conv1(x)
#         out = self.norm1(out)
#         out = self.relu1(out)

#         out = self.conv2(out)
#         out = self.norm2(out)

#         out = out + identity
#         out = self.relu2(out)
#         return out


# class UNet3DClassifier(nn.Module):
#     """
#     3D encoder classifier with residual blocks.
#     Input:  (B, 1, D, H, W)
#     Output: (B,) logits
#     """

#     def __init__(self, in_channels=1, base_channels=16, norm="group"):
#         super().__init__()

#         self.enc1 = ResidualBlock3D(in_channels, base_channels, norm=norm)
#         self.pool1 = nn.MaxPool3d(kernel_size=(1, 2, 2))   # keep depth at first stage

#         self.enc2 = ResidualBlock3D(base_channels, base_channels * 2, norm=norm)
#         self.pool2 = nn.MaxPool3d(kernel_size=(2, 2, 2))

#         self.enc3 = ResidualBlock3D(base_channels * 2, base_channels * 4, norm=norm)
#         self.pool3 = nn.MaxPool3d(kernel_size=(2, 2, 2))

#         self.bottleneck = ResidualBlock3D(base_channels * 4, base_channels * 8, norm=norm)

#         # self.global_avg_pool = nn.AdaptiveAvgPool3d(1)
#         # self.classifier = nn.Linear(base_channels * 8, 1)
#         self.global_avg_pool = nn.AdaptiveAvgPool3d(1)
#         self.global_max_pool = nn.AdaptiveMaxPool3d(1)
#         self.classifier = nn.Linear(base_channels * 16, 1)

#     def forward(self, x):
#         x = self.enc1(x)      # (B, C, 32, 512, 512)
#         x = self.pool1(x)     # (B, C, 32, 256, 256)

#         x = self.enc2(x)      # (B, 2C, 32, 256, 256)
#         x = self.pool2(x)     # (B, 2C, 16, 128, 128)

#         x = self.enc3(x)      # (B, 4C, 16, 128, 128)
#         x = self.pool3(x)     # (B, 4C, 8, 64, 64)

#         x = self.bottleneck(x)  # (B, 8C, 8, 64, 64)

#         # x = self.global_avg_pool(x)
#         # x = x.view(x.size(0), -1)
#         # x = self.classifier(x)
        
#         x_avg = self.global_avg_pool(x).view(x.size(0), -1)
#         x_max = self.global_max_pool(x).view(x.size(0), -1)
#         x = torch.cat([x_avg, x_max], dim=1)
#         x = self.classifier(x)
        
#         x = x.squeeze(1)

#         return x   

#%%

# class ResidualBlock3D(nn.Module):
#     def __init__(self, in_ch, out_ch, stride=1, groups=8):
#         super().__init__()

#         def gn(c):
#             g = min(groups, c)
#             while c % g != 0 and g > 1:
#                 g -= 1
#             return nn.GroupNorm(g, c)

#         self.conv1 = nn.Conv3d(in_ch, out_ch, kernel_size=3, stride=stride, padding=1, bias=False)
#         self.norm1 = gn(out_ch)
#         self.act1 = nn.SiLU(inplace=True)

#         self.conv2 = nn.Conv3d(out_ch, out_ch, kernel_size=3, padding=1, bias=False)
#         self.norm2 = gn(out_ch)

#         if in_ch != out_ch or stride != 1:
#             self.skip = nn.Sequential(
#                 nn.Conv3d(in_ch, out_ch, kernel_size=1, stride=stride, bias=False),
#                 gn(out_ch),
#             )
#         else:
#             self.skip = nn.Identity()

#         self.act2 = nn.SiLU(inplace=True)

#     def forward(self, x):
#         identity = self.skip(x)

#         out = self.conv1(x)
#         out = self.norm1(out)
#         out = self.act1(out)

#         out = self.conv2(out)
#         out = self.norm2(out)

#         out = out + identity
#         out = self.act2(out)
#         return out


# class ProjectionHead2D(nn.Module):
#     def __init__(self, in_ch, mid_ch=64):
#         super().__init__()
#         self.net = nn.Sequential(
#             nn.Conv2d(in_ch, mid_ch, kernel_size=3, padding=1, bias=False),
#             nn.GroupNorm(8 if mid_ch >= 8 else 1, mid_ch),
#             nn.SiLU(inplace=True),

#             nn.Conv2d(mid_ch, mid_ch, kernel_size=3, padding=1, bias=False),
#             nn.GroupNorm(8 if mid_ch >= 8 else 1, mid_ch),
#             nn.SiLU(inplace=True),
#         )
#         self.avg = nn.AdaptiveAvgPool2d(1)
#         self.max = nn.AdaptiveMaxPool2d(1)

#     def forward(self, x2d):
#         x2d = self.net(x2d)
#         a = self.avg(x2d).flatten(1)
#         m = self.max(x2d).flatten(1)
#         return torch.cat([a, m], dim=1)


# class WholeVolumeWakeNet(nn.Module):
#     def __init__(self, in_channels=2, base=16, dropout=0.25):
#         super().__init__()

#         self.stem = nn.Sequential(
#             nn.Conv3d(in_channels, base, kernel_size=3, padding=1, bias=False),
#             nn.GroupNorm(4 if base >= 4 else 1, base),
#             nn.SiLU(inplace=True),
#         )

#         self.block1 = ResidualBlock3D(base, base)
#         self.down1  = nn.Conv3d(base, base * 2, kernel_size=3, stride=(1, 2, 2), padding=1, bias=False)

#         self.block2 = ResidualBlock3D(base * 2, base * 2)
#         self.down2  = nn.Conv3d(base * 2, base * 4, kernel_size=3, stride=(2, 2, 2), padding=1, bias=False)

#         self.block3 = ResidualBlock3D(base * 4, base * 4)
#         self.down3  = nn.Conv3d(base * 4, base * 8, kernel_size=3, stride=(2, 2, 2), padding=1, bias=False)

#         self.block4 = ResidualBlock3D(base * 8, base * 8)
        
#         self.drop3d_1 = nn.Dropout3d(p=0.05)
#         self.drop3d_2 = nn.Dropout3d(p=0.10)

#         # 3D pooled branch
#         self.avg3d = nn.AdaptiveAvgPool3d(1)
#         self.max3d = nn.AdaptiveMaxPool3d(1)

#         # projection-aware 2D branch
#         self.proj_head = ProjectionHead2D(in_ch=base * 16, mid_ch=base * 4)

#         # final classifier
#         feat_3d = base * 8 * 2
#         feat_2d = base * 4 * 2
#         self.classifier = nn.Sequential(
#             nn.Linear(feat_3d + feat_2d, base * 8),
#             nn.SiLU(inplace=True),
#             nn.Dropout(dropout),
#             nn.Linear(base * 8, 1),
#         )
#         with torch.no_grad():
#             self.classifier[-1].bias.fill_(-1.0)

#     def forward(self, x):
#         x = self.stem(x)

#         x = self.block1(x)
#         x = self.down1(x)

#         x = self.block2(x)
#         x = self.down2(x)

#         # x = self.block3(x)
#         # x = self.down3(x)

#         # x = self.block4(x)
        
#         x = self.block3(x)
#         x = self.drop3d_1(x)
#         x = self.down3(x)
        
#         x = self.block4(x)
#         x = self.drop3d_2(x)

#         # 3D branch
#         f3a = self.avg3d(x).flatten(1)
#         f3m = self.max3d(x).flatten(1)
#         f3 = torch.cat([f3a, f3m], dim=1)

#         # projection-aware branch: use learned feature projections
#         x_mean = x.mean(dim=2)           # mean over depth
#         x_max  = x.max(dim=2).values     # max over depth
#         x2d = torch.cat([x_mean, x_max], dim=1)
#         f2 = self.proj_head(x2d)

#         out = self.classifier(torch.cat([f3, f2], dim=1)).squeeze(1)
#         return out 


class SubvolumeEncoder3D(nn.Module):
    """
    Shared encoder applied to each subvolume.
    Input:  (B*8, C, 4, H, W)
    Output: (B*8, feat_dim)
    """

    def __init__(self, in_channels=2, base=16, feat_dim=128):
        super().__init__()

        self.net = nn.Sequential(
            nn.Conv3d(in_channels, base, kernel_size=3, padding=1, bias=False),
            nn.GroupNorm(4 if base >= 4 else 1, base),
            nn.SiLU(inplace=True),

            nn.Conv3d(base, base * 2, kernel_size=3, stride=(1, 2, 2), padding=1, bias=False),
            nn.GroupNorm(8 if base * 2 >= 8 else 1, base * 2),
            nn.SiLU(inplace=True),

            nn.Conv3d(base * 2, base * 4, kernel_size=3, stride=(2, 2, 2), padding=1, bias=False),
            nn.GroupNorm(8 if base * 4 >= 8 else 1, base * 4),
            nn.SiLU(inplace=True),

            nn.Conv3d(base * 4, base * 4, kernel_size=3, padding=1, bias=False),
            nn.GroupNorm(8 if base * 4 >= 8 else 1, base * 4),
            nn.SiLU(inplace=True),
        )

        self.avg_pool = nn.AdaptiveAvgPool3d(1)
        self.max_pool = nn.AdaptiveMaxPool3d(1)

        self.proj = nn.Sequential(
            nn.Linear(base * 4 * 2, feat_dim),
            nn.SiLU(inplace=True),
        )

    def forward(self, x):
        x = self.net(x)
        x_avg = self.avg_pool(x).flatten(1)
        x_max = self.max_pool(x).flatten(1)
        x = torch.cat([x_avg, x_max], dim=1)
        x = self.proj(x)
        return x


class MILSubvolumeWakeNet(nn.Module):
    """
    Stage-1 MIL/attention model.

    Input:
        X: (B, C, 32, H, W)

    Output:
        vol_logit:     (B,)
        subvol_logits: (B, 8)
        attn_weights:  (B, 8)
    """

    def __init__(
        self,
        in_channels=2,
        n_subvolumes=8,
        subvol_depth=4,
        base=16,
        feat_dim=128,
        dropout=0.3,
    ):
        super().__init__()

        self.n_subvolumes = n_subvolumes
        self.subvol_depth = subvol_depth

        self.encoder = SubvolumeEncoder3D(
            in_channels=in_channels,
            base=base,
            feat_dim=feat_dim,
        )

        self.subvol_head = nn.Sequential(
            nn.Linear(feat_dim, 64),
            nn.SiLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

        self.attention = nn.Sequential(
            nn.Linear(feat_dim, 64),
            nn.Tanh(),
            nn.Linear(64, 1, bias=False),
        )

        self.vol_classifier = nn.Sequential(
            nn.Linear(feat_dim, 128),
            nn.SiLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(128, 1),
        )

        with torch.no_grad():
            self.subvol_head[-1].bias.fill_(-1.0)
            self.vol_classifier[-1].bias.fill_(-1.0)

    def forward(self, x, subvol_prior):
        B, C, D, H, W = x.shape

        expected_D = self.n_subvolumes * self.subvol_depth
        assert D == expected_D, f"Expected D={expected_D}, got D={D}"

        # (B, C, 32, H, W) -> (B, 8, C, 4, H, W)
        x = x.view(B, C, self.n_subvolumes, self.subvol_depth, H, W)
        x = x.permute(0, 2, 1, 3, 4, 5).contiguous()

        # (B, 8, C, 4, H, W) -> (B*8, C, 4, H, W)
        x = x.view(B * self.n_subvolumes, C, self.subvol_depth, H, W)

        # shared encoder
        feats = self.encoder(x)  # (B*8, feat_dim)

        # local subvolume logits
        subvol_logits = self.subvol_head(feats).view(B, self.n_subvolumes)

        # reshape features back to volume structure
        feats = feats.view(B, self.n_subvolumes, -1)  # (B, 8, feat_dim)

        # # attention over subvolumes
        # attn_logits = self.attention(feats).squeeze(-1)  # (B, 8)
        # attn_weights = torch.softmax(attn_logits, dim=1)
        attn_logits = self.attention(feats).squeeze(-1)          # (B, 8)
        prior_logits = torch.log(subvol_prior.clamp(min=1e-8))  # (B, 8)
        attn_weights = torch.softmax(attn_logits + prior_logits, dim=1)

        # weighted feature aggregation
        vol_feat = torch.sum(attn_weights.unsqueeze(-1) * feats, dim=1)

        # global volume logit
        vol_logit = self.vol_classifier(vol_feat).squeeze(1)

        return vol_logit, subvol_logits, attn_weights
    
    
    
#%%




# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# model = UNet3DClassifier(
#     in_channels=2,
#     base_channels=16,
#     norm="group"
# ).to(device)

# criterion = nn.BCEWithLogitsLoss()


# # # # this should be uncommented (if not on debug)

# # optimizer = optim.Adam(model.parameters(), lr=3e-4)

# # # this should be uncommented (if on debug)
# # optimizer = optim.Adam(model.parameters(), lr=3e-4)

# optimizer = optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)

# class AsymmetricFPLoss(nn.Module):
#     def __init__(self, gamma_pos=0.0, gamma_neg=3.0, clip=0.05, eps=1e-8):
#         super().__init__()
#         self.gamma_pos = gamma_pos
#         self.gamma_neg = gamma_neg
#         self.clip = clip
#         self.eps = eps

#     def forward(self, logits, targets):
#         targets = targets.float()
#         probs = torch.sigmoid(logits)

#         pos_probs = probs
#         neg_probs = 1.0 - probs

#         if self.clip is not None and self.clip > 0:
#             neg_probs = torch.clamp(neg_probs + self.clip, max=1.0)

#         pos_loss = targets * torch.log(pos_probs.clamp(min=self.eps))
#         neg_loss = (1.0 - targets) * torch.log(neg_probs.clamp(min=self.eps))

#         if self.gamma_pos > 0 or self.gamma_neg > 0:
#             pos_weight = torch.pow(1.0 - pos_probs, self.gamma_pos) * targets
#             neg_weight = torch.pow(1.0 - neg_probs, self.gamma_neg) * (1.0 - targets)
#             loss = pos_weight * pos_loss + neg_weight * neg_loss
#         else:
#             loss = pos_loss + neg_loss

#         return -loss.mean()
    
    
    

#%%
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# model = WholeVolumeWakeNet(in_channels=2, base=16, dropout=0.5).to(device)
# optimizer = torch.optim.AdamW(
#     model.parameters(),
#     lr=1e-4,
#     weight_decay=1e-4
# )


model = MILSubvolumeWakeNet(
    in_channels=2,
    n_subvolumes=8,
    subvol_depth=4,
    base=16,
    feat_dim=128,
    dropout=0.3,
).to(device)



optimizer = torch.optim.AdamW(
    model.parameters(),
    lr=1e-4,
    weight_decay=1e-4
)

# # criterion = nn.BCEWithLogitsLoss()
criterion = nn.BCEWithLogitsLoss(
    pos_weight=torch.tensor([POS_WEIGHT], device=device)
)
lambda_sub = 0.3


# criterion = AsymmetricFPLoss(
#     gamma_pos=0.0,
#     gamma_neg=3.0,
#     clip=0.05
# )



#%%
model.eval()

# for X, y in train_dataloader:
#     X = X.to(device)
#     y = y.to(device).float()

#     with torch.no_grad():
#         logits = model(X)
#         probs = torch.sigmoid(logits)

#     print("\n--- Sanity check ---")
#     print("X shape:", X.shape)
#     print("y shape:", y.shape, y.dtype)
#     print("logits shape:", logits.shape)
#     print("logits:", logits)
#     print("probs:", probs)
#     break

for X, vol_label, subvol_label in train_dataloader:
    X = X.to(device)
    vol_label = vol_label.to(device).float()
    subvol_label = subvol_label.to(device).float()

    with torch.no_grad():
        # vol_logit, subvol_logits, attn_weights = model(X)
        subvol_prior = subvol_label.clone()
        subvol_prior = subvol_prior.clamp(min=0)
        subvol_prior = subvol_prior / (subvol_prior.sum(dim=1, keepdim=True) + 1e-8)
        subvol_prior = subvol_prior.to(device)
        
        vol_logit, subvol_logits, attn_weights = model(X, subvol_prior)
        probs = torch.sigmoid(vol_logit)

    print("\n--- Sanity check ---")
    print("X shape:", X.shape)
    print("vol_label shape:", vol_label.shape, vol_label.dtype)
    print("subvol_label shape:", subvol_label.shape, subvol_label.dtype)
    print("vol_logit shape:", vol_logit.shape)
    print("subvol_logits shape:", subvol_logits.shape)
    print("attn_weights shape:", attn_weights.shape)
    print("vol probs:", probs)
    print("attention weights:", attn_weights)

    break




#%%

scheduler = torch.optim.lr_scheduler.OneCycleLR(
    optimizer,
    max_lr=2e-4,
    epochs=num_epochs,
    steps_per_epoch=len(train_dataloader),
    pct_start=0.15,
    div_factor=10,
    final_div_factor=20,
)


#%%


# def train_one_epoch(model, dataloader, criterion, optimizer, device, scheduler):
#     model.train()
#     running_loss = 0.0
#     running_correct = 0
#     total = 0
    
#     grad_norm_sum = 0.0
#     grad_steps = 0

#     for X, y in dataloader:
#         X = X.to(device)
#         y = y.to(device).float()

#         optimizer.zero_grad()

#         logits = model(X)
#         loss = criterion(logits, y)

#         loss.backward()
        
#         total_grad_sq = 0.0
#         for p in model.parameters():
#             if p.grad is not None:
#                 g = p.grad.detach()
#                 total_grad_sq += g.norm(2).item() ** 2
        
#         grad_norm = total_grad_sq ** 0.5
#         grad_norm_sum += grad_norm
#         grad_steps += 1
        
#         optimizer.step()
#         scheduler.step()

#         probs = torch.sigmoid(logits)
#         preds = (probs > 0.5).float()

#         running_loss += loss.item() * X.size(0)
#         running_correct += (preds == y).sum().item()
#         total += X.size(0)

#     epoch_loss = running_loss / total if total > 0 else float("nan")
#     epoch_acc = running_correct / total if total > 0 else float("nan")
#     mean_grad_norm = grad_norm_sum / grad_steps if grad_steps > 0 else float("nan")

#     return epoch_loss, epoch_acc, mean_grad_norm

# def train_one_epoch(model, dataloader, criterion, optimizer, device, accum_steps=4):
#     model.train()
#     running_loss = 0.0
#     running_correct = 0
#     total = 0

#     grad_norm_sum = 0.0
#     grad_steps = 0

#     optimizer.zero_grad()

#     for step, (X, y) in enumerate(dataloader):
#         X = X.to(device)
#         y = y.to(device).float()

#         # ---- train-time augmentation: flips only ----
#         # X shape is [B, C, D, H, W]
#         if torch.rand(1).item() < 0.5:
#             X = torch.flip(X, dims=[3])   # flip H
#         if torch.rand(1).item() < 0.5:
#             X = torch.flip(X, dims=[4])   # flip W

#         logits = model(X)
#         loss = criterion(logits, y)

#         # keep full loss value for logging
#         running_loss += loss.item() * X.size(0)

#         # gradient accumulation
#         loss = loss / accum_steps
#         loss.backward()

#         # gradient norm diagnostic
#         total_grad_sq = 0.0
#         for p in model.parameters():
#             if p.grad is not None:
#                 g = p.grad.detach()
#                 total_grad_sq += g.norm(2).item() ** 2
#         grad_norm = total_grad_sq ** 0.5
#         grad_norm_sum += grad_norm
#         grad_steps += 1
        
#         torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)


#         if ((step + 1) % accum_steps == 0) or ((step + 1) == len(dataloader)):
#             optimizer.step()
#             optimizer.zero_grad()

#         probs = torch.sigmoid(logits)
#         preds = (probs > 0.5).float()

#         running_correct += (preds == y).sum().item()
#         total += X.size(0)

#     epoch_loss = running_loss / total if total > 0 else float("nan")
#     epoch_acc = running_correct / total if total > 0 else float("nan")
#     mean_grad_norm = grad_norm_sum / grad_steps if grad_steps > 0 else float("nan")

#     return epoch_loss, epoch_acc, mean_grad_norm

def train_one_epoch(
    model,
    dataloader,
    criterion,
    optimizer,
    device,
    lambda_sub=0.3,
    accum_steps=4,
):
    model.train()

    running_loss = 0.0
    running_vol_loss = 0.0
    running_sub_loss = 0.0
    running_correct = 0
    total = 0

    grad_norm_sum = 0.0
    grad_steps = 0

    optimizer.zero_grad()

    for step, (X, vol_label, subvol_label) in enumerate(dataloader):
        X = X.to(device)
        vol_label = vol_label.to(device).float()
        subvol_label = subvol_label.to(device).float()

        # train-time augmentation: flips only
        if torch.rand(1).item() < 0.5:
            X = torch.flip(X, dims=[3])
        if torch.rand(1).item() < 0.5:
            X = torch.flip(X, dims=[4])

        # vol_logit, subvol_logits, attn_weights = model(X)
        # Build prior from subvol_label before calling model
        subvol_prior = subvol_label.clone().clamp(min=0)
        subvol_prior = subvol_prior / (subvol_prior.sum(dim=1, keepdim=True) + 1e-8)        
        vol_logit, subvol_logits, attn_weights = model(X, subvol_prior)

        loss_vol = criterion(vol_logit, vol_label)
        loss_sub = criterion(subvol_logits, subvol_label)

        loss = loss_vol + lambda_sub * loss_sub

        running_loss += loss.item() * X.size(0)
        running_vol_loss += loss_vol.item() * X.size(0)
        running_sub_loss += loss_sub.item() * X.size(0)

        loss = loss / accum_steps
        loss.backward()

        total_grad_sq = 0.0
        for p in model.parameters():
            if p.grad is not None:
                g = p.grad.detach()
                total_grad_sq += g.norm(2).item() ** 2

        grad_norm = total_grad_sq ** 0.5
        grad_norm_sum += grad_norm
        grad_steps += 1

        if ((step + 1) % accum_steps == 0) or ((step + 1) == len(dataloader)):
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            optimizer.zero_grad()
            scheduler.step()

        probs = torch.sigmoid(vol_logit)
        preds = (probs > 0.5).float()

        running_correct += (preds == vol_label).sum().item()
        total += X.size(0)

    epoch_loss = running_loss / total if total > 0 else float("nan")
    epoch_vol_loss = running_vol_loss / total if total > 0 else float("nan")
    epoch_sub_loss = running_sub_loss / total if total > 0 else float("nan")
    epoch_acc = running_correct / total if total > 0 else float("nan")
    mean_grad_norm = grad_norm_sum / grad_steps if grad_steps > 0 else float("nan")

    return epoch_loss, epoch_vol_loss, epoch_sub_loss, epoch_acc, mean_grad_norm





# def evaluate(model, dataloader, criterion, device):
#     model.eval()
#     running_loss = 0.0
#     running_correct = 0
#     total = 0

#     with torch.no_grad():
#         for X, y in dataloader:
#             X = X.to(device)
#             y = y.to(device).float()

#             logits = model(X)
#             loss = criterion(logits, y)

#             probs = torch.sigmoid(logits)
#             preds = (probs > 0.5).float()

#             running_loss += loss.item() * X.size(0)
#             running_correct += (preds == y).sum().item()
#             total += X.size(0)

#     epoch_loss = running_loss / total if total > 0 else float("nan")
#     epoch_acc = running_correct / total if total > 0 else float("nan")
#     return epoch_loss, epoch_acc

# def evaluate(model, dataloader, criterion, device, return_details=False):
#     model.eval()
#     running_loss = 0.0
#     running_correct = 0
#     total = 0

#     all_probs = []
#     all_logits = []
#     all_y = []

#     with torch.no_grad():
#         for X, y in dataloader:
#             X = X.to(device)
#             y = y.to(device).float()

#             logits = model(X)
#             loss = criterion(logits, y)

#             probs = torch.sigmoid(logits)
#             preds = (probs > 0.5).float()

#             running_loss += loss.item() * X.size(0)
#             running_correct += (preds == y).sum().item()
#             total += X.size(0)

#             all_probs.extend(probs.detach().cpu().numpy().tolist())
#             all_logits.extend(logits.detach().cpu().numpy().tolist())
#             all_y.extend(y.detach().cpu().numpy().tolist())

#     epoch_loss = running_loss / total if total > 0 else float("nan")
#     epoch_acc = running_correct / total if total > 0 else float("nan")

#     if return_details:
#         return (
#             epoch_loss,
#             epoch_acc,
#             np.array(all_probs),
#             np.array(all_logits),
#             np.array(all_y),
#         )

#     return epoch_loss, epoch_acc


def evaluate(
    model,
    dataloader,
    criterion,
    device,
    lambda_sub=0.3,
    return_details=False,
):
    model.eval()

    running_loss = 0.0
    running_vol_loss = 0.0
    running_sub_loss = 0.0
    running_correct = 0
    total = 0

    all_probs = []
    all_logits = []
    all_y = []
    all_attn = []

    with torch.no_grad():
        for X, vol_label, subvol_label in dataloader:
            X = X.to(device)
            vol_label = vol_label.to(device).float()
            subvol_label = subvol_label.to(device).float()

            # vol_logit, subvol_logits, attn_weights = model(X)
            subvol_prior = subvol_label.clone().clamp(min=0)
            subvol_prior = subvol_prior / (subvol_prior.sum(dim=1, keepdim=True) + 1e-8)
            
            vol_logit, subvol_logits, attn_weights = model(X, subvol_prior)

            loss_vol = criterion(vol_logit, vol_label)
            loss_sub = criterion(subvol_logits, subvol_label)
            loss = loss_vol + lambda_sub * loss_sub

            probs = torch.sigmoid(vol_logit)
            preds = (probs > 0.5).float()

            running_loss += loss.item() * X.size(0)
            running_vol_loss += loss_vol.item() * X.size(0)
            running_sub_loss += loss_sub.item() * X.size(0)
            running_correct += (preds == vol_label).sum().item()
            total += X.size(0)

            all_probs.extend(probs.detach().cpu().numpy().tolist())
            all_logits.extend(vol_logit.detach().cpu().numpy().tolist())
            all_y.extend(vol_label.detach().cpu().numpy().tolist())
            all_attn.extend(attn_weights.detach().cpu().numpy().tolist())

    epoch_loss = running_loss / total if total > 0 else float("nan")
    epoch_vol_loss = running_vol_loss / total if total > 0 else float("nan")
    epoch_sub_loss = running_sub_loss / total if total > 0 else float("nan")
    epoch_acc = running_correct / total if total > 0 else float("nan")

    if return_details:
        return (
            epoch_loss,
            epoch_vol_loss,
            epoch_sub_loss,
            epoch_acc,
            np.array(all_probs),
            np.array(all_logits),
            np.array(all_y),
            np.array(all_attn),
        )

    return epoch_loss, epoch_vol_loss, epoch_sub_loss, epoch_acc


def threshold_stats(y_true, probs, thr):
    y_true = np.asarray(y_true).astype(int)
    probs = np.asarray(probs)

    y_pred = (probs > thr).astype(int)

    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))

    recall = tp / (tp + fn + 1e-12)
    precision = tp / (tp + fp + 1e-12)
    acc = (tp + tn) / (tp + tn + fp + fn + 1e-12)
    fn_over_tp = fn / (tp + 1e-12)

    return {
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "recall": float(recall),
        "precision": float(precision),
        "acc": float(acc),
        "fn_over_tp": float(fn_over_tp),
    }


def find_best_threshold_for_precision_priority(
    y_true,
    probs,
    thresholds=THRESH_GRID,
    min_recall=0.05
):
    best_thr = DEFAULT_THRESHOLD
    best_stats = None

    for thr in thresholds:
        stats = threshold_stats(y_true, probs, thr)

        # avoid trivial all-negative solution
        if stats["recall"] < min_recall:
            continue

        if best_stats is None:
            best_thr = thr
            best_stats = stats
            continue

        better = False

        # primary goal: minimize FP/TP
        if stats["fp"] / (stats["tp"] + 1e-12) < best_stats["fp"] / (best_stats["tp"] + 1e-12):
            better = True

        # tie-break 1: higher precision
        elif np.isclose(
            stats["fp"] / (stats["tp"] + 1e-12),
            best_stats["fp"] / (best_stats["tp"] + 1e-12)
        ):
            if stats["precision"] > best_stats["precision"]:
                better = True

            # tie-break 2: more TP
            elif np.isclose(stats["precision"], best_stats["precision"]):
                if stats["tp"] > best_stats["tp"]:
                    better = True

        if better:
            best_thr = thr
            best_stats = stats

    # fallback if no threshold satisfies min_recall
    if best_stats is None:
        best_thr = DEFAULT_THRESHOLD
        best_stats = threshold_stats(y_true, probs, best_thr)

    return best_thr, best_stats

#%%
# num_epochs = 1



train_losses, train_accs = [], []
val_losses, val_accs = [], []

best_val_loss = float("inf")

# for epoch in range(num_epochs):
#     train_loss, train_acc = train_one_epoch(model, train_dataloader, criterion, optimizer, device)
#     val_loss, val_acc = evaluate(model, valid_dataloader, criterion, device)

#     train_losses.append(train_loss)
#     train_accs.append(train_acc)
#     val_losses.append(val_loss)
#     val_accs.append(val_acc)

#     print(f"Epoch {epoch+1}/{num_epochs}")
#     print(f"  Train loss: {train_loss:.4f} | Train acc: {train_acc:.4f}")
#     print(f"  Val   loss: {val_loss:.4f} | Val   acc: {val_acc:.4f}")

best_fp_over_tp = float("inf")
best_precision = -1.0
best_tp = -1
best_epoch = -1
best_threshold = DEFAULT_THRESHOLD
best_path = "best_model_precision.pt"



for epoch in range(num_epochs):
    # train_loss, train_acc, mean_grad_norm = train_one_epoch(
    # model, train_dataloader, criterion, optimizer, device, accum_steps=4
    # )
    
    # val_loss, val_acc, val_probs, val_logits, val_y = evaluate(
    #     model, valid_dataloader, criterion, device, return_details=True
    # )
    train_loss, train_vol_loss, train_sub_loss, train_acc, mean_grad_norm = train_one_epoch(
    model,
    train_dataloader,
    criterion,
    optimizer,
    device,
    lambda_sub=lambda_sub,
    accum_steps=4,
    )
    
    val_loss, val_vol_loss, val_sub_loss, val_acc, val_probs, val_logits, val_y, val_attn = evaluate(
        model,
        valid_dataloader,
        criterion,
        device,
        lambda_sub=lambda_sub,
        return_details=True,
    )

    try:
        val_auc = roc_auc_score(val_y, val_probs)
    except Exception:
        val_auc = float("nan")
        
    # improved = False

    # if val_loss < best_val_loss:
    #     best_val_loss = val_loss
    #     improved = True

    # if val_acc > best_val_acc:
    #     best_val_acc = val_acc
    #     best_epoch = epoch + 1
    #     improved = True

    # if improved:
    #     torch.save(model.state_dict(), best_path)
    #     epochs_no_improve = 0
    # else:
    #     epochs_no_improve += 1
        
        # choose threshold that prioritizes precision 
    thr_epoch, thr_stats = find_best_threshold_for_precision_priority(
    val_y,
    val_probs,
    thresholds=THRESH_GRID,
    min_recall=0.20
)

    current_fp_over_tp = thr_stats["fp"] / (thr_stats["tp"] + 1e-12)
    current_precision = thr_stats["precision"]
    current_tp = thr_stats["tp"]
    
    improved = False
    
    if best_epoch == -1:
        improved = True
    elif current_fp_over_tp < best_fp_over_tp:
        improved = True
    elif np.isclose(current_fp_over_tp, best_fp_over_tp):
        if current_precision > best_precision:
            improved = True
        elif np.isclose(current_precision, best_precision):
            if current_tp > best_tp:
                improved = True
    
    if improved:
        best_fp_over_tp = current_fp_over_tp
        best_precision = current_precision
        best_tp = current_tp
        best_epoch = epoch + 1
        best_threshold = thr_epoch
        torch.save(model.state_dict(), best_path)

        

    train_losses.append(train_loss)
    train_accs.append(train_acc)
    val_losses.append(val_loss)
    val_accs.append(val_acc)

    current_fp_over_tp = thr_stats["fp"] / (thr_stats["tp"] + 1e-12)

    print(f"Epoch {epoch+1}/{num_epochs}")
    # print(f"  Train loss: {train_loss:.4f} | Train acc: {train_acc:.4f}")
    # print(f"  Val   loss: {val_loss:.4f} | Val   acc@0.5: {val_acc:.4f}")
    print(f"  Train total loss: {train_loss:.4f} | vol: {train_vol_loss:.4f} | subvol: {train_sub_loss:.4f} | acc: {train_acc:.4f}")
    print(f"  Val   total loss: {val_loss:.4f} | vol: {val_vol_loss:.4f} | subvol: {val_sub_loss:.4f} | acc@0.5: {val_acc:.4f}")
    print("  Mean attention weights:", np.mean(val_attn, axis=0))
    print(f"  Val AUC: {val_auc:.4f}")
    print(f"  Mean grad norm: {mean_grad_norm:.6f}")
    print(f"  Best precision-priority threshold this epoch: {thr_epoch:.2f}")
    print(
        f"  TP={thr_stats['tp']} FN={thr_stats['fn']} "
        f"FP={thr_stats['fp']} TN={thr_stats['tn']}"
    )
    print(
        f"  Recall={thr_stats['recall']:.4f} | "
        f"Precision={thr_stats['precision']:.4f} | "
        f"FP/TP={current_fp_over_tp:.4f}"
    )
    
    
#%%
current_fp_over_tp = thr_stats["fp"] / (thr_stats["tp"] + 1e-12)

print(f"Epoch {epoch+1}/{num_epochs}")
print(f"  Train loss: {train_loss:.4f} | Train acc: {train_acc:.4f}")
print(f"  Val   loss: {val_loss:.4f} | Val   acc@0.5: {val_acc:.4f}")
print(f"  Val AUC: {val_auc:.4f}")
print(f"  Mean grad norm: {mean_grad_norm:.6f}")
print(f"  Best precision-priority threshold this epoch: {thr_epoch:.2f}")
print(
    f"  TP={thr_stats['tp']} FN={thr_stats['fn']} "
    f"FP={thr_stats['fp']} TN={thr_stats['tn']}"
)
print(
    f"  Recall={thr_stats['recall']:.4f} | "
    f"Precision={thr_stats['precision']:.4f} | "
    f"FP/TP={current_fp_over_tp:.4f}"
)

