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
parser.add_argument('--slices_signal', type=int, default=32, help='')
# parser.add_argument('--percentage_positiveWakeSig', type=float, default=10, help='')
parser.add_argument('--validation_fraction', type=float, default=0.1, help='')
parser.add_argument('--train_tt_fraction', type=float, default=0.8, help='')
parser.add_argument('--wake_top_percentage', type=float, default=20, help='')
parser.add_argument('--void_percentage', type=float, default=0, help='')


parser.add_argument('--batch_size', type=int, default=32, help='')
parser.add_argument('--num_workers', type=int, default=0, help='')
parser.add_argument('--num_epochs', type=int, default=10, help='')

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
from torch.utils.data import DataLoader, Dataset, Subset
from torchvision import transforms
# import torch.utils.data as data

# from efficientnet_pytorch import EfficientNet
import torchvision.models as models

# Training function def
# import copy 

# data load

import os
from PIL import Image

# set the random seeds

import random

# math

import numpy as np


# Import the spliter function from DataCleaning.py'
import sys

# from DataCleaning import file_list,list_all_files,balanced_list_of_files
path_analy = os.getcwd() +'/' 
# path_analy = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path_analy)
from Data_preSelection  import  *

# keep track of true positives, false positives, true negatives, and false negatives for each class
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedShuffleSplit


# To count the number of images in each label
from collections import Counter


# Reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)


#%%

# Define the data transformations:

    
pretrained_means = [0.485, 0.456, 0.406]
pretrained_stds = [0.229, 0.224, 0.225]

train_transforms = transforms.Compose([
                           transforms.Resize(pretrained_size),
                           # transforms.RandomRotation(5),
                           transforms.RandomHorizontalFlip(0.5),
                           # transforms.RandomCrop(pretrained_size, padding=10),
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

# # Define the Custom Dataset Class

# class CustomImageDataset(Dataset):
#     def __init__(self, file_list, transform=None):
#         """
#         Args:
#             file_list (list): List of file paths to be included in the dataset.
#             transform (callable, optional): Optional transform to be applied on a sample.
#         """
#         self.file_list = file_list
#         self.transform = transform

#         # Extract classes from the parent of the immediate parent directory
#         self.classes = list(set([os.path.basename(os.path.dirname(os.path.dirname(file_path))) for file_path in file_list]))
#         self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}
#         self.labels = [self.class_to_idx[os.path.basename(os.path.dirname(os.path.dirname(file_path)))] for file_path in file_list]

#     def __len__(self):
#         return len(self.file_list)

#     def __getitem__(self, idx):
#         img_path = self.file_list[idx]
#         image = Image.open(img_path).convert("RGB")

#         # Extract class name from the parent of the immediate parent directory
#         class_name = os.path.basename(os.path.dirname(os.path.dirname(img_path)))
#         label = self.class_to_idx[class_name]

#         if self.transform:
#             image = self.transform(image)

#         return image, label

# class WakeBinDataset(Dataset):
#     """
#     Dataset that uses the *paths* in file_list only as metadata to locate
#     the corresponding .bin file (sample, anglid, sliceid, wake/no-wake),
#     and then reads the actual image slice from the .bin volume.

#     Assumes that each .bin has shape (n_slices, H, W) in float32.
#     """

#     def __init__(self, file_list, transform=None):
#         """
#         Args:
#             file_list (list[str]): List of PNG-like paths (only for metadata).
#             path_bin (str): Root directory where .bin files live.
#             n_slices (int): Number of slices stored in each .bin file.
#             img_size (int): Expected (H, W) size of each slice (square).
#             transform (callable, optional): Torchvision transforms.
#         """
#         self.file_list = file_list
#         # self.path_bin = path_bin
#         # self.n_slices = n_slices
#         # self.img_size = img_size
#         self.transform = transform

#         # Same class logic as before (wake / nowake folder name)
#         self.classes = sorted(
#             set(os.path.basename(os.path.dirname(os.path.dirname(fp)))
#                 for fp in file_list)
#         )
#         self.class_to_idx = {cls_name: idx for idx, cls_name in enumerate(self.classes)}
#         self.labels = [
#             self.class_to_idx[
#                 os.path.basename(os.path.dirname(os.path.dirname(fp)))
#             ]
#             for fp in file_list
#         ]

#         # Regex to decode sample / anglid / tilesize / sliceid from the png-like path
#         self.pattern = re.compile(
#             r"sample(\d+)-anglid_(\d+)-2dproj_z3_ts(\d+)_sl(\d+)\.png"
#         )

#     def __len__(self):
#         return len(self.file_list)

#     def _parse_meta(self, path):
#         """
#         Extract (sample, anglid, tilesize, sliceid, wake_info) from a file path.
#         """
#         m = self.pattern.search(path)
#         if not m:
#             raise ValueError(f"Could not parse metadata from path: {path}")
#         sample, anglid, tilesize, sliceid = map(int, m.groups())
#         wake_info = os.path.basename(os.path.dirname(os.path.dirname(path)))
#         return sample, anglid, tilesize, sliceid, wake_info

#     def __getitem__(self, idx):
#         meta_path = self.file_list[idx]
#         sample, anglid, tilesize, sliceid, wake_info = self._parse_meta(meta_path)

#         # Build .bin path. Adjust this to YOUR real layout.
#         # Example: <path_bin>/<wake_info>/sample<sample>/sample<sample>_anglid_<anglid>.bin
#         bin_dir = os.path.join(self.path_bin, wake_info, f"sample{sample}")
#         bin_path = os.path.join(bin_dir, f"sample{sample}_anglid_{anglid}.bin")

#         # Load the 3D volume from .bin
#         # Assumed: float32, shape = (n_slices, img_size, img_size)
#         vol = np.fromfile(bin_path, dtype=np.float32)

#         # You MUST adjust this reshape to match how you actually store the bin:
#         vol = vol.reshape(self.n_slices, self.img_size, self.img_size)

#         # Select the slice (sliceid is 1-based in your naming)
#         slice_idx = sliceid - 1
#         if not (0 <= slice_idx < self.n_slices):
#             raise IndexError(
#                 f"sliceid {sliceid} out of range for n_slices={self.n_slices} "
#                 f"(path: {meta_path})"
#             )
#         img2d = vol[slice_idx, :, :]  # shape (H, W)

#         # Normalize to [0,1] for safety; you can change this depending on your stats
#         img_min, img_max = img2d.min(), img2d.max()
#         if img_max > img_min:
#             img_norm = (img2d - img_min) / (img_max - img_min)
#         else:
#             img_norm = np.zeros_like(img2d, dtype=np.float32)

#         # Convert to uint8 grayscale 0–255, then to RGB PIL
#         img_uint8 = (img_norm * 255).astype(np.uint8)
#         img_pil = Image.fromarray(img_uint8, mode="L").convert("RGB")

#         if self.transform:
#             img_pil = self.transform(img_pil)

#         label = self.labels[idx]
#         return img_pil, label

#%%

# Define the Custom Dataset Class


# def load_bin_slice(Nmesh,BoxSize_,filename,redshift,slice_idx):

#     Nmesh_x, Nmesh_y, Nmesh_z = Nmesh
    
#     # Calculate the number of elements in a 2D slice (Nmesh_x * Nmesh_y)
#     elements_per_slice = Nmesh_x * Nmesh_y
    
#     # Calculate the offset in bytes for the desired slice
#     # Since each float32 is 4 bytes, multiply by 4 to get the byte offset
#     offset = slice_idx * elements_per_slice * 4  # offset in bytes
        
#     # Open the file and read only the necessary slice
#     with open(filename, 'rb') as file:
#         # Move the file pointer to the calculated offset
#         file.seek(offset)
        
#         # Read only the required number of elements for one slice
#         img2d = np.fromfile(file, dtype=np.float32, count=elements_per_slice)
#         # Reshape to 2D array (Nmesh_x, Nmesh_y)
        

#     return img2d



# def load_bin_slice(
#     filename: str,
#     grid_shape: tuple[int, int, int],
#     slice_id: int,
#     dtype: np.dtype = np.float32,
# ) -> np.ndarray:
#     """
#     Load a single 2D slice from a 3D volume stored in a raw binary (.bin) file.

#     Assumes the file contains a 3D array with shape:
#         (Nx, Ny, Nz)
#     stored in C order as:
#         slice 0 (all x,y), then slice 1, ..., slice Nz-1.

#     Parameters
#     ----------
#     filename : str
#         Path to the .bin file.
#     grid_shape : (int, int, int)
#         Tuple (Nx, Ny, Nz) giving the mesh size in each dimension.
#         We assume slices are planes of shape (Nx, Ny), indexed along z.
#     slice_id : int
#         1-based index of the slice along z to load (1 <= slice_id <= Nz).
#         This matches the “slice 1..slices_signal” convention.
#     dtype : np.dtype, optional
#         Data type stored in the file. Default is np.float32.

#     Returns
#     -------
#     img2d : np.ndarray
#         2D array of shape (Nx, Ny) containing the requested slice.
#     """

#     Nx, Ny, Nz = grid_shape

#     # --- Validate slice_id (1-based) ---
#     if not (1 <= slice_id <= Nz):
#         raise ValueError(
#             f"slice_id={slice_id} out of range; expected 1..{Nz}"
#         )

#     # Convert to 0-based index internally
#     slice_idx = slice_id - 1

#     # Number of elements in one 2D slice
#     elements_per_slice = Nx * Ny

#     # Size of each element in bytes
#     itemsize = np.dtype(dtype).itemsize

#     # Offset in bytes from the beginning of the file
#     offset_bytes = slice_idx * elements_per_slice * itemsize

#     # Expected total file size (for sanity check)
#     expected_size_bytes = elements_per_slice * Nz * itemsize

#     with open(filename, "rb") as f:
#         # Go to end to check file size
#         f.seek(0, 2)
#         file_size = f.tell()
#         if file_size != expected_size_bytes:
#             raise ValueError(
#                 f"File {filename} has size {file_size} bytes, "
#                 f"expected {expected_size_bytes} for grid_shape={grid_shape} "
#                 f"and dtype={dtype}."
#             )

#         # Seek to the beginning of the desired slice
#         f.seek(offset_bytes, 0)

#         # Read exactly one slice worth of data
#         flat = np.fromfile(f, dtype=dtype, count=elements_per_slice)

#     if flat.size != elements_per_slice:
#         raise ValueError(
#             f"Could not read full slice from {filename}: "
#             f"got {flat.size} elements, expected {elements_per_slice}."
#         )

#     # Reshape to (Nx, Ny)
#     img2d = flat.reshape(Nx, Ny)
#     img2d = img2d.T

#     return img2d


def compute_file_means(
    bin_paths: list[str],
    grid_shape: tuple[int, int, int],
    dtype: np.dtype = np.float32,
) -> dict[str, float]:
    """
    Lê cada arquivo .bin da lista **uma única vez** e calcula
    a média global do volume 3D (Nx, Ny, Nz).

    Retorna:
        { caminho_bin : média_global }
    """
    Nx, Ny, Nz = grid_shape
    n_elems = Nx * Ny * Nz
    file_means = {}

    # usa set() para não recalcular a média de arquivos repetidos
    unique_paths = sorted(set(bin_paths))

    for path in unique_paths:
        flat = np.fromfile(path, dtype=dtype, count=n_elems)
        if flat.size != n_elems:
            raise ValueError(
                f"Arquivo {path} tem {flat.size} elementos, "
                f"mas esperados {n_elems} para grid_shape={grid_shape}."
            )
        file_means[path] = float(flat.mean())

    return file_means



def load_bin_slice(
    filename: str,
    grid_shape: tuple[int, int, int],
    slice_id: int,
    mean_global: float | None = None,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    """
    Lê um único slice 2D de um volume 3D armazenado em arquivo .bin bruto e,
    opcionalmente, aplica a transformação estilo Matlab:
        dc = (map - mean_global) / mean_global
        map_3d_slices = atan((dc + 1) * 16)

    filename   : caminho do .bin
    grid_shape : (Nx, Ny, Nz)
    slice_id   : índice 1-based ao longo de z (1..Nz)
    mean_global: média do volume correspondente (já pré-computada)
    dtype      : tipo dos dados no arquivo (float32, etc)
    """
    Nx, Ny, Nz = grid_shape

    if not (1 <= slice_id <= Nz):
        raise ValueError(f"slice_id={slice_id} fora do intervalo 1..{Nz}")

    slice_idx = slice_id - 1
    elements_per_slice = Nx * Ny
    itemsize = np.dtype(dtype).itemsize
    offset_bytes = slice_idx * elements_per_slice * itemsize
    expected_size_bytes = elements_per_slice * Nz * itemsize

    with open(filename, "rb") as f:
        f.seek(0, 2)
        file_size = f.tell()
        if file_size != expected_size_bytes:
            raise ValueError(
                f"Arquivo {filename} tem {file_size} bytes, "
                f"mas esperados {expected_size_bytes} para grid_shape={grid_shape}, dtype={dtype}"
            )
        f.seek(offset_bytes, 0)
        flat = np.fromfile(f, dtype=dtype, count=elements_per_slice)

    if flat.size != elements_per_slice:
        raise ValueError(
            f"Não foi possível ler o slice completo de {filename}: "
            f"lidos {flat.size}, esperados {elements_per_slice} elementos."
        )

    img2d = flat.reshape(Nx, Ny).T  # mesmo transpose que você já estava usando

    # --- transformação estilo Matlab, se média for fornecida ---
    if mean_global is not None and mean_global != 0.0:
        # dc = (map - mean_global) / mean_global
        dc = (img2d - mean_global) / mean_global
        # map_3d_slices = atan((dc + 1) * 16)
        img2d = np.arctan((dc + 1.0) * 16.0).astype(np.float32)

    return img2d



## verify if it loads correctly


class BinSliceDataset(Dataset):
    """
    Dataset para slices armazenados em volumes .bin (Nx, Ny, Nz).

    Cada amostra é definida por:
      - bin_paths[i]   : caminho do arquivo .bin
      - slice_ids[i]   : índice 1-based do slice dentro do volume
      - wake_infos[i]  : string indicando wake / nowake

    O dataset:
      - lê o slice 2D correto dentro do arquivo .bin
      - aplica a transformação tipo Matlab (usando média pré-computada do volume)
      - normaliza para [0,1], converte para PIL RGB
      - devolve (imagem_transformada, label)

    Atributos
    ---------
    labels : list[int]
        Labels inteiros (0 = nowake, 1 = wake), um por amostra.
    classes : list[str]
        Lista de nomes de classe, indexed por label.
    class_to_idx : dict[str, int]
        Mapeia nomes de classe para índices inteiros.
    """

    def __init__(
        self,
        bin_paths: list[str],
        slice_ids: list[int],
        wake_infos: list[str],
        grid_shape: tuple[int, int, int],
        file_means: dict[str, float],
        transform=None,
        dtype: np.dtype = np.float32,
    ):
        if not (len(bin_paths) == len(slice_ids) == len(wake_infos)):
            raise ValueError(
                "bin_paths, slice_ids e wake_infos devem ter o mesmo comprimento."
            )

        self.bin_paths  = bin_paths
        self.slice_ids  = slice_ids
        self.wake_infos = wake_infos
        self.grid_shape = tuple(grid_shape)
        self.transform  = transform
        self.dtype      = dtype

        # dicionário { caminho_bin : média_global_do_volume }
        self.file_means = file_means

        # Mesmas constantes que você já usa
        self.class_to_idx = {
            NOWAKE_NAME: 0,
            WAKE_NAME: 1,
        }

        # classes[idx] -> nome
        self.classes = [None] * len(self.class_to_idx)
        for name, idx in self.class_to_idx.items():
            self.classes[idx] = name

        # labels inteiros pré-computados
        self.labels = [self._label_from_wake_info(wi) for wi in wake_infos]

    def __len__(self):
        return len(self.bin_paths)

    def _label_from_wake_info(self, wi: str) -> int:
        if WAKE_NAME in wi:
            return 1
        elif NOWAKE_NAME in wi:
            return 0
        else:
            raise ValueError(f"wake_info desconhecido: {wi}")

    def __getitem__(self, idx):
        bin_path  = self.bin_paths[idx]
        slice_id  = self.slice_ids[idx]   # 1-based
        label     = self.labels[idx]

        mean_global = self.file_means.get(bin_path, None)
        if mean_global is None:
            raise KeyError(f"Média não encontrada em file_means para arquivo: {bin_path}")

        # --- lê slice 2D pronto (já com atan((dc+1)*16)) ---
        img2d = load_bin_slice(
            filename=bin_path,
            grid_shape=self.grid_shape,
            slice_id=slice_id,
            mean_global=mean_global,
            dtype=self.dtype,
        )

        # --- Normalização [0,1] por slice (igual ao caso png) ---
        mn, mx = img2d.min(), img2d.max()
        if mx > mn:
            img_norm = (img2d - mn) / (mx - mn)
        else:
            img_norm = np.zeros_like(img2d, dtype=np.float32)

        img_uint8 = (img_norm * 255).astype(np.uint8)
        img_pil = Image.fromarray(img_uint8, mode="L").convert("RGB")

        if self.transform is not None:
            img_out = self.transform(img_pil)
        else:
            img_out = img_pil

        return img_out, label

#%%

# Pipeline summary:
# 1. Load precomputed void/signal statistics
# 2. List all .bin files and compute one global mean per file
# 3. Expand each 3D file into one entry per slice
# 4. Split by unique sample to avoid leakage between validation and train/test
# 5. Select extreme wake/no-wake cases using signal and void criteria
# 6. Build PyTorch datasets and dataloaders


#%%


# Load datasets


all_data_nowake_void, all_data_wake_void = data_void_out(rang, n_angle, slices_void, wake_spec, path_void)
all_data_nowake_signal, all_data_wake_signal = data_signal_out(rang, n_angle, slices_signal, wake_spec, path_WakeSignal)


files_list_all = list_all_files(path_data,rang)
file_means = compute_file_means(files_list_all, grid_shape=Nmesh, dtype=np.float32)

samples_all, anglids_all, wake_infos_all = extract_info_bin(files_list_all)


files_slices_all    = []
samples_slices_all    = []
anglids_slices_all    = []
wake_infos_slices_all = []
sliceids_all          = []

for file, sample, anglid, wake_info in zip(files_list_all,samples_all, anglids_all, wake_infos_all):
    for slice_id in range(1, slices_signal + 1):
        files_slices_all.append(file)
        samples_slices_all.append(sample)
        anglids_slices_all.append(anglid)
        wake_infos_slices_all.append(wake_info)
        sliceids_all.append(slice_id)
        
        
# range_start = rang[0]

data_void = extract_stat(samples_slices_all, anglids_slices_all, sliceids_all, wake_infos_slices_all, all_data_nowake_void, all_data_wake_void, rang)
data_signal = extract_stat(samples_slices_all, anglids_slices_all, sliceids_all, wake_infos_slices_all, all_data_nowake_signal, all_data_wake_signal, rang)
data_signal_diff = extract_stat_diff(samples_slices_all, anglids_slices_all, sliceids_all, wake_infos_slices_all, all_data_nowake_signal, all_data_wake_signal, rang)


# split val and trainTest

# first split per sample (so the validation is independent)

validation_indices, train_test_indices = split_unique_samples(samples_slices_all,wake_infos_slices_all, validation_fraction)

files_list_validation = [files_slices_all[i] for i in validation_indices]
files_list_trainTest = [files_slices_all[i] for i in train_test_indices]

samples_val    = [samples_slices_all[i] for i in validation_indices]
anglids_val    = [anglids_slices_all[i] for i in validation_indices]
sliceids_val   = [sliceids_all[i] for i in validation_indices]
wake_infos_val = [wake_infos_slices_all[i] for i in validation_indices]


samples_tt    = [samples_slices_all[i] for i in train_test_indices]
anglids_tt    = [anglids_slices_all[i] for i in train_test_indices]
sliceids_tt   = [sliceids_all[i] for i in train_test_indices]
wake_infos_tt = [wake_infos_slices_all[i] for i in train_test_indices]

#%%


# range_start = rang[0]
data_void_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_void, all_data_wake_void, rang)
data_signal_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, rang)
data_signal_diff_val = extract_stat_diff(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, rang)

# range_start = rang[0]
data_void_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_void, all_data_wake_void, rang)
data_signal_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, rang)
data_signal_diff_tt = extract_stat_diff(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, rang)



# # Select top wake values, keep void_percentage = 50%
# selected_files, selected_positions, selected_signal_diff = select_extreme_files(
#     data_signal_diff, wake_infos_all, files_list_all, data_void,
#     wake_top_percentage, void_percentage
# )




# validation

# Select top 50% wake values, keep void_percentage 
selected_files_val, selected_positions_val, selected_signal_diff_val = select_extreme_files(
    data_signal_diff_val, wake_infos_val, files_list_validation, data_void_val,
    wake_top_percentage, void_percentage
)

find_extreme_files(selected_signal_diff_val, selected_files_val)
selected_sliceids_val   = [sliceids_val[i] for i in selected_positions_val]
selected_wake_infos_val = [wake_infos_val[i] for i in selected_positions_val]




# train and test

# Select top 50% wake values, keep void_percentage 
selected_files_tt, selected_positions_tt, selected_signal_diff_tt = select_extreme_files(
    data_signal_diff_tt, wake_infos_tt, files_list_trainTest, data_void_tt,
    wake_top_percentage, void_percentage
)

selected_sliceids_tt   = [sliceids_tt[i] for i in selected_positions_tt]
selected_wake_infos_tt = [wake_infos_tt[i] for i in selected_positions_tt]


find_extreme_files(selected_signal_diff_tt, selected_files_tt)


#%%

# shape of your 3D grid in each bin file:
# grid_shape = (slices_signal, original_H, original_W)
# e.g. if you know it's 32x512x512:
# grid_shape = (slices_signal, 512, 512)

valid_data__ = BinSliceDataset(
    bin_paths=selected_files_val,
    slice_ids=selected_sliceids_val,
    wake_infos=selected_wake_infos_val,
    grid_shape=Nmesh,
    file_means=file_means,
    transform=test_transforms,
)

test_train_data__ = BinSliceDataset(
    bin_paths=selected_files_tt,
    slice_ids=selected_sliceids_tt,
    wake_infos=selected_wake_infos_tt,
    grid_shape=Nmesh,
    file_means=file_means,
    transform=test_transforms,  # or train_transforms if you prefer
)


# Full validation set: ALL entries belonging to the validation samples
# (before select_extreme_files filtering)
valid_all_data__ = BinSliceDataset(
    bin_paths=files_list_validation,
    slice_ids=sliceids_val,
    wake_infos=wake_infos_val,
    grid_shape=Nmesh,
    file_means=file_means,
    transform=test_transforms,
)




#%%

# valid_data__ = CustomImageDataset(file_list=selected_files_val, transform=test_transforms)
# test_train_data__ = CustomImageDataset(file_list=selected_files_tt, transform=test_transforms)



n_train_examples = int(len(test_train_data__) * train_tt_fraction)
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

# train_data_, test_data_ = data.random_split(test_train_data__,
#                                             [n_train_examples, n_test_examples])



#%%




# Create data loaders.

valid_dataloader = DataLoader(valid_data__, batch_size=batch_size, shuffle=False, num_workers=num_workers)
test_dataloader = DataLoader(test_data_, batch_size=batch_size, shuffle=False, num_workers=num_workers)
train_dataloader = DataLoader(train_data_, batch_size=batch_size, shuffle=True, num_workers=num_workers)

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




valid_all_dataloader = DataLoader(
    valid_all_data__,
    batch_size=batch_size,
    shuffle=False,
    num_workers=num_workers
)





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
    
# Unfreeze last feature block
for p in model.features[-1].parameters():
    p.requires_grad = True    
    
# Define the loss function and optimizer
criterion = nn.BCEWithLogitsLoss()  # Binary Cross Entropy with Logits Loss
optimizer = optim.Adam(
    [
        {"params": model.classifier[-1].parameters()},
        {"params": model.classifier[-2].parameters()},
        {"params": model.features[-1].parameters()},
    ],
    lr=1e-4,   # ↓ a bit to avoid blowing things up
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


for inputs, labels in train_dataloader:
    print("DEBUG first batch:", inputs.shape, inputs.dtype, inputs.device)
    break
#%%

# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# model = model.to(device)


# import torch.nn.functional as F

def prepare_binary_logits(outputs: torch.Tensor) -> torch.Tensor:
    """
    Convert model output to shape [batch] for BCEWithLogitsLoss.
    Handles the special case batch_size == 1.
    """
    if outputs.shape == torch.Size([1, 1]):
        return outputs.squeeze().unsqueeze(0)
    return outputs.squeeze()



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
            # # outputs = outputs.squeeze()
            # if outputs.shape == torch.Size([1, 1]):  # outputs is a scalar           
            #     outputs = outputs.squeeze().unsqueeze(0)   # Reshape scalar to [1]                
            # else:
            #     outputs = outputs.squeeze()    # Squeeze the output if necessary
            outputs = prepare_binary_logits(outputs)

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
        cm = confusion_matrix(all_labels, all_preds, labels=[0, 1])
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
        cm = confusion_matrix(all_labels, all_preds, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        valid_class_0_accuracy = tn / (tn + fp)
        valid_class_1_accuracy = tp / (tp + fn)

        print(f'Epoch {epoch + 1}/{num_epochs} - Validation Loss: {valid_epoch_loss:.4f} Acc: {valid_epoch_acc:.4f}')
        print(f'Validation Class 0 Accuracy: {valid_class_0_accuracy:.4f}, Validation Class 1 Accuracy: {valid_class_1_accuracy:.4f}')

    return model


def evaluate_full_validation(
    model,
    dataloader,
    sample_ids,
    wake_infos,
    criterion=None,
    dataset_name="ALL validation entries"
):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.eval()
    model.to(device)

    all_labels = []
    all_preds = []
    all_probs = []
    running_loss = 0.0

    with torch.no_grad():
        for inputs, labels in dataloader:
            inputs = inputs.to(device)
            labels = labels.to(device).float()

            outputs = model(inputs)

            if outputs.shape == torch.Size([1, 1]):
                outputs = outputs.squeeze().unsqueeze(0)
            else:
                outputs = outputs.squeeze()

            probs = torch.sigmoid(outputs)
            preds = (probs > 0.5).long()

            if criterion is not None:
                loss = criterion(outputs, labels)
                running_loss += loss.item() * inputs.size(0)

            all_labels.append(labels.cpu().numpy().astype(int))
            all_preds.append(preds.cpu().numpy().astype(int))
            all_probs.append(probs.cpu().numpy())

    all_labels = np.concatenate(all_labels)
    all_preds = np.concatenate(all_preds)
    all_probs = np.concatenate(all_probs)

    if len(sample_ids) != len(all_labels):
        raise ValueError(
            f"sample_ids has length {len(sample_ids)}, but predictions have length {len(all_labels)}"
        )
    if len(wake_infos) != len(all_labels):
        raise ValueError(
            f"wake_infos has length {len(wake_infos)}, but predictions have length {len(all_labels)}"
        )

    # ---------- Slice-level metrics ----------
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()

    total = cm.sum()
    acc = (tn + tp) / total if total > 0 else float("nan")
    class_0_acc = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    class_1_acc = tp / (tp + fn) if (tp + fn) > 0 else float("nan")

    print(f"\n===== {dataset_name}: slice-level results =====")
    if criterion is not None:
        avg_loss = running_loss / len(dataloader.dataset)
        print(f"Loss: {avg_loss:.4f}")
    print(f"Accuracy: {acc:.4f}")
    print(f"Class 0 Accuracy: {class_0_acc:.4f}")
    print(f"Class 1 Accuracy: {class_1_acc:.4f}")
    print("Confusion matrix [ [TN, FP], [FN, TP] ]:")
    print(cm)

    # ---------- Group by (sample_id, wake_info) ----------
    grouped = {}
    for sid, wi, y_true, prob, y_pred in zip(sample_ids, wake_infos, all_labels, all_probs, all_preds):
        key = (sid, wi)
        if key not in grouped:
            grouped[key] = {
                "labels": [],
                "probs": [],
                "preds": [],
            }
        grouped[key]["labels"].append(int(y_true))
        grouped[key]["probs"].append(float(prob))
        grouped[key]["preds"].append(int(y_pred))

    print(f"\n===== {dataset_name}: grouped by (sample_id, wake_info) =====")
    grouped_true = []
    grouped_pred = []

    for (sid, wi), vals in sorted(grouped.items()):
        mean_prob = float(np.mean(vals["probs"]))
        frac_positive = float(np.mean(vals["preds"]))
        true_label = int(round(np.mean(vals["labels"])))
        pred_from_mean = int(mean_prob > 0.5)
        n_entries = len(vals["labels"])

        grouped_true.append(true_label)
        grouped_pred.append(pred_from_mean)

        print(
            f"sample={sid} | wake_info={wi} | n_entries={n_entries} | "
            f"true={true_label} | mean_prob={mean_prob:.4f} | "
            f"frac_pred_positive={frac_positive:.4f} | pred_from_mean={pred_from_mean}"
        )

    grouped_true = np.array(grouped_true)
    grouped_pred = np.array(grouped_pred)
    grouped_acc = np.mean(grouped_true == grouped_pred) if len(grouped_true) > 0 else float("nan")

    print(f"\nGrouped accuracy over (sample_id, wake_info): {grouped_acc:.4f}")

    return {
        "slice_labels": all_labels,
        "slice_preds": all_preds,
        "slice_probs": all_probs,
        "grouped": grouped,
        "slice_confusion_matrix": cm,
        "slice_accuracy": acc,
        "grouped_accuracy": grouped_acc,
    }


#%%


# Call the training function
trained_model = train_model(model, train_dataloader, valid_dataloader, criterion, optimizer, num_epochs=num_epochs)


#%%

full_val_results = evaluate_full_validation(
    trained_model,
    valid_all_dataloader,
    sample_ids=samples_val,
    wake_infos=wake_infos_val,
    criterion=criterion,
    dataset_name="ALL entries from validation samples"
)
# for inputs, labels in train_dataloader:
#     print(inputs)
#     print(labels)
