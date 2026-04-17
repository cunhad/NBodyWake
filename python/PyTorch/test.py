#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:52:37 2026

@author: asus
"""

#%%



## compare bin and png


import numpy as np



# same function used in the bin nn analisis

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


import numpy as np

def compute_volume_mean(filename, grid_shape, dtype=np.float32):
    Nx, Ny, Nz = grid_shape
    total_elems = Nx * Ny * Nz
    data = np.fromfile(filename, dtype=dtype, count=total_elems)
    if data.size != total_elems:
        raise ValueError("Bad file size")
    return data.mean()



def load_bin_slice(
    filename: str,
    grid_shape: tuple[int, int, int],
    slice_id: int,
    mean_global: float,
    dtype: np.dtype = np.float32,
) -> np.ndarray:
    Nx, Ny, Nz = grid_shape
    if not (1 <= slice_id <= Nz):
        raise ValueError(...)
    slice_idx = slice_id - 1
    elements_per_slice = Nx * Ny
    itemsize = np.dtype(dtype).itemsize
    offset_bytes = slice_idx * elements_per_slice * itemsize

    with open(filename, "rb") as f:
        f.seek(offset_bytes, 0)
        flat = np.fromfile(f, dtype=dtype, count=elements_per_slice)

    if flat.size != elements_per_slice:
        raise ValueError("Could not read full slice")

    img2d = flat.reshape(Nx, Ny).T  # same orientation as before

    # Apply MATLAB-style transform using precomputed mean
    dc = (img2d - mean_global) / mean_global
    img2d_out = np.arctan((dc + 1.0) * 16.0)
    return img2d_out





import matplotlib.pyplot as plt
from PIL import Image

png_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs/4Mpc_2048c_1024p_zi63_nowakem/sample5001/sample5001-anglid_1-2dproj_z3_ts32_sl1.png"
bin_path = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_nowakem/sample5001/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin"
grid_shape = (512, 512, 32)
slice_id   = 1   # same “sl1” as in the PNG

# 1) Load PNG
img_png = Image.open(png_path).convert("L")
img_png = np.array(img_png, dtype=np.float32)
#resize stuff
img_png = Image.fromarray(img_png)
img_png = img_png.resize((512, 512), Image.BILINEAR)


# 2) Load BIN slice with your function
mean_global = compute_volume_mean(bin_path, grid_shape)
img_bin = load_bin_slice(bin_path, grid_shape, slice_id, mean_global)
# img_bin = (img_bin - img_bin.min()) / (img_bin.max() - img_bin.min() + 1e-8)

fig, ax = plt.subplots(1, 2, figsize=(8,4))
ax[0].imshow(img_png, cmap="gray")
ax[0].set_title("PNG")
ax[1].imshow(img_bin, cmap="gray")
ax[1].set_title("BIN slice")
plt.show()