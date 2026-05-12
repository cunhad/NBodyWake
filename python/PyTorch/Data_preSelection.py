#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:51:39 2025

@author: asus
"""




#%%

# Void Data



# def data_void_out(rang, n_angle, slices, wake_spec, path):
    
#     import numpy as np


#     # slices_dp = 8
    
#     # Initialize an empty list to store the flattened data
#     all_data_nowake = np.empty((len(rang), n_angle, slices))
#     all_data_wake = np.empty((len(rang), n_angle, slices))
    
    
    
#     # for simul in range(3001,3010):
#     # for simul in range(5001,5100):
#     # for simul in range(3001,3001+1):
#     for i, simul in enumerate(rang):    
#         # print("sample"+str(simul))
#         # filename_nowake = path_void+wake_spec[0]+"sample"+str(simul)+"_2ds4t3_curv_z3_stat.txt" 
#         filename_nowake = path + wake_spec[0] + f"sample{simul}_2d_void_z3_stat.txt"
#         filename_wake =   path+wake_spec[1]+"sample"+str(simul)+"_2d_void_z3_stat.txt" 
        
#         try:
#             data_array_nowake = np.loadtxt(filename_nowake, delimiter='\t')
#             all_data_nowake[i, :, :] = data_array_nowake
#             # Check if the array has the expected dimensions
#             if data_array_nowake.shape != (n_angle, slices):
#                 print(f"The array no wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_nowake.shape}.")
#                 break
#             # Flatten the data array to 1D and append to the list
#             # all_data_nowake.extend(data_array_nowake.flatten())
            
#         except Exception as e:
#             print(f"Error loading file {filename_nowake}: {e}")
#             break
        
#         try:
#             data_array_wake = np.loadtxt(filename_wake, delimiter='\t')
#             all_data_wake[i, :, :] = data_array_wake
#             # Check if the array has the expected dimensions
#             if data_array_wake.shape != (n_angle, slices):
#                 print(f"The array wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_wake.shape}.")
#                 break
#             # Flatten the data array to 1D and append to the list
#             # all_data_wake.extend(data_array_wake.flatten())
            
#         except Exception as e:
#             print(f"Error loading file {filename_nowake}: {e}")
#             break
        
#     return all_data_nowake, all_data_wake
    
def data_void_out(rang, n_angle, slices, wake_spec, path):
    import numpy as np
    import os

    all_data_nowake = np.full((len(rang), n_angle, slices), np.nan, dtype=float)
    all_data_wake   = np.full((len(rang), n_angle, slices), np.nan, dtype=float)

    def load_void_file(base_path):
        """
        Try standard filename first, then down1 variant.
        """
        if os.path.isfile(base_path):
            return np.loadtxt(base_path, delimiter="\t")

        alt_path = base_path.replace("_z3_stat.txt", "_z3_down1_stat.txt")
        if os.path.isfile(alt_path):
            return np.loadtxt(alt_path, delimiter="\t")

        raise FileNotFoundError(
            f"Neither file exists:\n  {base_path}\n  {alt_path}"
        )

    for i, simul in enumerate(rang):

        filename_nowake = (
            path + wake_spec[0] + f"sample{simul}_2d_void_z3_stat.txt"
        )
        filename_wake = (
            path + wake_spec[1] + f"sample{simul}_2d_void_z3_stat.txt"
        )

        try:
            data_array_nowake = load_void_file(filename_nowake)
            if data_array_nowake.shape != (n_angle, slices):
                raise ValueError(
                    f"No-wake sample {simul} has shape {data_array_nowake.shape}, "
                    f"expected {(n_angle, slices)}"
                )
            all_data_nowake[i, :, :] = data_array_nowake

        except Exception as e:
            print(f"Error loading no-wake file for sample {simul}: {e}")
            continue

        try:
            data_array_wake = load_void_file(filename_wake)
            if data_array_wake.shape != (n_angle, slices):
                raise ValueError(
                    f"Wake sample {simul} has shape {data_array_wake.shape}, "
                    f"expected {(n_angle, slices)}"
                )
            all_data_wake[i, :, :] = data_array_wake

        except Exception as e:
            print(f"Error loading wake file for sample {simul}: {e}")
            continue

    return all_data_nowake, all_data_wake



# def data_void_out3d(rang, n_angle, slices, wake_spec, path):
#     import numpy as np
#     import os

#     all_data_nowake = np.full((len(rang), n_angle), np.nan, dtype=float)
#     all_data_wake   = np.full((len(rang), n_angle), np.nan, dtype=float)

#     def load_void_file(base_path):
#         """
#         Try standard filename first, then down1 variant.
#         """
#         if os.path.isfile(base_path):
#             return np.loadtxt(base_path, delimiter="\t")

#         alt_path = base_path.replace("_z3_stat.txt", "_z3_down1_stat.txt")
#         if os.path.isfile(alt_path):
#             return np.loadtxt(alt_path, delimiter="\t")

#         raise FileNotFoundError(
#             f"Neither file exists:\n  {base_path}\n  {alt_path}"
#         )

#     for i, simul in enumerate(rang):

#         filename_nowake = (
#             path + wake_spec[0] + f"sample{simul}_2d_void_z3_stat.txt"
#         )
#         filename_wake = (
#             path + wake_spec[1] + f"sample{simul}_2d_void_z3_stat.txt"
#         )

#         try:
#             data_array_nowake = load_void_file(filename_nowake)
#             if data_array_nowake.shape != (n_angle, slices):
#                 raise ValueError(
#                     f"No-wake sample {simul} has shape {data_array_nowake.shape}, "
#                     f"expected {(n_angle, slices)}"
#                 )
#             all_data_nowake[i, :] = data_array_nowake[:,slices-1]

#         except Exception as e:
#             print(f"Error loading no-wake file for sample {simul}: {e}")
#             break

#         try:
#             data_array_wake = load_void_file(filename_wake)
#             if data_array_wake.shape != (n_angle, slices):
#                 raise ValueError(
#                     f"Wake sample {simul} has shape {data_array_wake.shape}, "
#                     f"expected {(n_angle, slices)}"
#                 )
#             all_data_wake[i, :] = data_array_wake[:,slices-1]

#         except Exception as e:
#             print(f"Error loading wake file for sample {simul}: {e}")
#             break

#     return all_data_nowake, all_data_wake

def data_void_out3d(rang, n_angle, slices, wake_spec, path):
    import numpy as np
    import os

    # Fill with NaN so failed loads do not leave garbage values
    all_data_nowake = np.full((len(rang), n_angle), np.nan, dtype=float)
    all_data_wake   = np.full((len(rang), n_angle), np.nan, dtype=float)

    def load_void_file(base_path):
        """
        Try standard filename first, then down1 variant.
        """
        if os.path.isfile(base_path):
            return np.loadtxt(base_path, delimiter="\t")

        alt_path = base_path.replace("_z3_stat.txt", "_z3_down1_stat.txt")
        if os.path.isfile(alt_path):
            return np.loadtxt(alt_path, delimiter="\t")

        raise FileNotFoundError(
            f"Neither file exists:\n  {base_path}\n  {alt_path}"
        )

    for i, simul in enumerate(rang):

        filename_nowake = (
            path + wake_spec[0] + f"sample{simul}_2d_void_z3_stat.txt"
        )
        filename_wake = (
            path + wake_spec[1] + f"sample{simul}_2d_void_z3_stat.txt"
        )

        # ---------------- no-wake ----------------
        try:
            data_array_nowake = load_void_file(filename_nowake)

            if data_array_nowake.shape != (n_angle, slices):
                raise ValueError(
                    f"No-wake sample {simul} has shape {data_array_nowake.shape}, "
                    f"expected {(n_angle, slices)}"
                )

            # Keep the same convention you were using before:
            # take the last slice statistic for each angle
            all_data_nowake[i, :] = data_array_nowake[:, slices - 1]

        except Exception as e:
            print(f"Error loading no-wake file for sample {simul}: {e}")
            # Leave this row as NaN and continue

        # ---------------- wake ----------------
        try:
            data_array_wake = load_void_file(filename_wake)

            if data_array_wake.shape != (n_angle, slices):
                raise ValueError(
                    f"Wake sample {simul} has shape {data_array_wake.shape}, "
                    f"expected {(n_angle, slices)}"
                )

            all_data_wake[i, :] = data_array_wake[:, slices - 1]

        except Exception as e:
            print(f"Error loading wake file for sample {simul}: {e}")
            # Leave this row as NaN and continue

    print("NaNs in all_data_nowake_void:", np.isnan(all_data_nowake).sum())
    print("NaNs in all_data_wake_void:", np.isnan(all_data_wake).sum())

    return all_data_nowake, all_data_wake

#%%

# Signal Data
def data_signal_out(rang, n_angle, slices, wake_spec, path):
    
    import numpy as np


        
    
    # Initialize an empty list to store the flattened data
    all_data_nowake = np.full((len(rang), n_angle, slices), np.nan, dtype=float)
    all_data_wake   = np.full((len(rang), n_angle, slices), np.nan, dtype=float)
    
    
    
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
                continue
            # Flatten the data array to 1D and append to the list
            # all_data_nowake.extend(data_array_nowake.flatten())
            
        except Exception as e:
            print(f"Error loading file {filename_nowake}: {e}")
            continue
        
        try:
            data_array_wake = np.loadtxt(filename_wake, delimiter='\t')
            all_data_wake[i, :, :] = data_array_wake
            # Check if the array has the expected dimensions
            if data_array_wake.shape != (n_angle, slices):
                print(f"The array wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_wake.shape}.")
                continue
            # Flatten the data array to 1D and append to the list
            # all_data_wake.extend(data_array_wake.flatten())
            
        except Exception as e:
            print(f"Error loading file {filename_nowake}: {e}")
            continue
        
    return all_data_nowake, all_data_wake

# def data_signal_out3d(rang, n_angle, slices, wake_spec, path):
    
#     import numpy as np


        
    
#     # Initialize an empty list to store the flattened data
#     all_data_nowake = np.full((len(rang), n_angle), np.nan, dtype=float)
#     all_data_wake   = np.full((len(rang), n_angle), np.nan, dtype=float)
    
    
    
#     # for simul in range(3001,3010):
#     # for simul in range(5001,5100):
#     # for simul in range(3001,3001+1):
#     for i, simul in enumerate(rang):    
#         # print("sample"+str(simul))
#         # filename_nowake = path+wake_spec[0]+"sample"+str(simul)+"_2ds4t3_curv_z3_stat.txt" 
#         filename_nowake = path + wake_spec[0] + f"sample{simul}_2ds4t3dp_curv_z3_stat.txt"
#         filename_wake =   path+wake_spec[1]+"sample"+str(simul)+"_2ds4t3dp_curv_z3_stat.txt" 
        
#         try:
#             data_array_nowake = np.loadtxt(filename_nowake, delimiter='\t')
            
#             # Check if the array has the expected dimensions
#             if data_array_nowake.shape != (n_angle, slices):
#                 print(f"The array no wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_nowake.shape}.")
#                 break
            
#             # Sum over slices → shape (n_angle,)
#             all_data_nowake[i, :] = data_array_nowake.sum(axis=1)
#             # Flatten the data array to 1D and append to the list
#             # all_data_nowake.extend(data_array_nowake.flatten())
            
#         except Exception as e:
#             print(f"Error loading file {filename_nowake}: {e}")
#             break
        
#         try:
#             data_array_wake = np.loadtxt(filename_wake, delimiter='\t')
            
#             # Check if the array has the expected dimensions
#             if data_array_wake.shape != (n_angle, slices):
#                 print(f"The array wake {simul} does not have {n_angle} rows and {slices} columns. Its shape is {data_array_wake.shape}.")
#                 break
#             all_data_wake[i, :] = data_array_wake.sum(axis=1)
#             # Flatten the data array to 1D and append to the list
#             # all_data_wake.extend(data_array_wake.flatten())
            
#         except Exception as e:
#             print(f"Error loading file {filename_nowake}: {e}")
#             break
        
#     return all_data_nowake, all_data_wake


def data_signal_out3d(rang, n_angle, slices, wake_spec, path):
    import numpy as np

    all_data_nowake = np.full((len(rang), n_angle), np.nan, dtype=float)
    all_data_wake   = np.full((len(rang), n_angle), np.nan, dtype=float)

    for i, simul in enumerate(rang):
        filename_nowake = path + wake_spec[0] + f"sample{simul}_2ds4t3dp_curv_z3_stat.txt"
        filename_wake   = path + wake_spec[1] + f"sample{simul}_2ds4t3dp_curv_z3_stat.txt"

        try:
            data_array_nowake = np.loadtxt(filename_nowake, delimiter="\t")
            if data_array_nowake.shape == (n_angle, slices):
                all_data_nowake[i, :] = data_array_nowake.sum(axis=1)
            else:
                print(f"Bad no-wake shape for sample {simul}: {data_array_nowake.shape}")
        except Exception as e:
            print(f"Error loading file {filename_nowake}: {e}")

        try:
            data_array_wake = np.loadtxt(filename_wake, delimiter="\t")
            if data_array_wake.shape == (n_angle, slices):
                all_data_wake[i, :] = data_array_wake.sum(axis=1)
            else:
                print(f"Bad wake shape for sample {simul}: {data_array_wake.shape}")
        except Exception as e:
            print(f"Error loading file {filename_wake}: {e}")

    return all_data_nowake, all_data_wake

def data_signal_subvol_out3d(rang, n_angle, slices, wake_spec, path):
    import numpy as np

    all_data_nowake = np.full((len(rang), n_angle, slices), np.nan, dtype=float)
    all_data_wake   = np.full((len(rang), n_angle, slices), np.nan, dtype=float)

    for i, simul in enumerate(rang):
        filename_nowake = path + wake_spec[0] + f"sample{simul}_2ds4t3dp_curv_z3_stat.txt"
        filename_wake   = path + wake_spec[1] + f"sample{simul}_2ds4t3dp_curv_z3_stat.txt"

        try:
            arr_nowake = np.loadtxt(filename_nowake, delimiter="\t")
            if arr_nowake.shape == (n_angle, slices):
                all_data_nowake[i, :, :] = arr_nowake
            else:
                print(f"Bad no-wake subvol shape for sample {simul}: {arr_nowake.shape}")
        except Exception as e:
            print(f"Error loading subvol no-wake file {filename_nowake}: {e}")

        try:
            arr_wake = np.loadtxt(filename_wake, delimiter="\t")
            if arr_wake.shape == (n_angle, slices):
                all_data_wake[i, :, :] = arr_wake
            else:
                print(f"Bad wake subvol shape for sample {simul}: {arr_wake.shape}")
        except Exception as e:
            print(f"Error loading subvol wake file {filename_wake}: {e}")

    return all_data_nowake, all_data_wake

#%%

# Other functions


# def list_all_files(folder_path):
    
#     import os
    
#     all_files = []
    
#     # os.walk generates the file names in a directory tree, including all subdirectories
#     for root, dirs, files in os.walk(folder_path):
#         # 🔥 PRUNE directories: only descend into folders containing "wake"
#         # dirs[:] = [d for d in dirs if "wake" in d]
#         for file in files:
#             # Construct full file path
#             file_path = os.path.join(root, file)
#             all_files.append(file_path)
    
#     return all_files

# def list_all_files(root_path):
#     import os

#     all_files = []

#     for d in os.listdir(root_path):
#         if "wake" not in d:
#             continue

#         full_dir = os.path.join(root_path, d)
#         if not os.path.isdir(full_dir):
#             continue

#         for root, _, files in os.walk(full_dir):
#             for file in files:
#                 all_files.append(os.path.join(root, file))

#     return all_files

def list_all_files(root_path, samples=None):
    """
    List all files under root_path, restricted to:
      - directories whose name contains 'wake' (e.g. nowakem, wakeGmu...)
      - optionally, only files belonging to samples in `samples`
        (we detect sample ID via 'sampleXXXX' in the path).

    Parameters
    ----------
    root_path : str
        Top-level directory (PNG or bin root).
    samples : list[int] or None
        If provided, only keep files whose path contains 'sampleXXXX'
        with XXXX in this list. Use e.g. your `rang`.

    Returns
    -------
    all_files : list[str]
        Full file paths.
    """
    import os
    import re

    all_files = []
    sample_set = set(samples) if samples is not None else None

    for d in os.listdir(root_path):
        full_dir = os.path.join(root_path, d)
        if not os.path.isdir(full_dir):
            continue

        # Keep your original behaviour: only enter sim dirs that contain "wake"
        # (this matches both 'nowakem' and 'wakeGmu4t10m8zi10m', and skips 'plots_9', etc.)
        if "wake" not in d:
            continue

        for root, _, files in os.walk(full_dir):
            for file in files:
                full_path = os.path.join(root, file)

                # If no sample filter, behave like your original version
                if sample_set is None:
                    all_files.append(full_path)
                    continue

                # Otherwise, require 'sampleXXXX' and check if XXXX in samples
                m = re.search(r"sample(\d+)", full_path)
                if not m:
                    continue

                sample_id = int(m.group(1))
                if sample_id in sample_set:
                    all_files.append(full_path)

    return all_files


#%%




def extract_info(paths: list[str]) -> tuple[list[int], list[int], list[int], list[int], list[str]]:
    
    import re
    import os
    
    
    pattern = r"sample(\d+)-anglid_(\d+)-2dproj_z3_ts(\d+)_sl(\d+)\.png"
    
    samples, anglids, tilesizes, sliceids, wake_infos = [], [], [], [], []
    
    for path in paths:
        match = re.search(pattern, path)
        if match:
            sample, anglid, tilesize, sliceid = map(int, match.groups())
            samples.append(sample)
            anglids.append(anglid)
            tilesizes.append(tilesize)
            sliceids.append(sliceid)
            
            # Extract folder name
            folder_name = os.path.basename(os.path.dirname(os.path.dirname(path)))
            wake_infos.append(folder_name)
    
    return samples, anglids, tilesizes, sliceids, wake_infos



def extract_info_bin(paths: list[str]) -> tuple[list[int], list[int], list[str]]:
    """
    Extract sample id, angular id, and wake info from .bin file paths.

    Assumes paths like (no-wake):
      .../4Mpc_2048c_1024p_zi63_nowakem/sample5001/data/.../anglid_1/.../*.bin

    and like (wake, with extra folder under sample):
      .../4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/half_lin_cutoff_.../data/.../anglid_2/.../*.bin

    Returns:
      samples   : list[int]
      anglids   : list[int]
      wake_infos: list[str]
    """

    import re
    import os

    samples = []
    anglids = []
    wake_infos = []

    sample_pattern = re.compile(r"sample(\d+)")
    anglid_pattern = re.compile(r"anglid_(\d+)")

    for path in paths:
        norm_path = os.path.normpath(path)
        parts = norm_path.split(os.sep)

        # --- sample id ---
        m_sample = sample_pattern.search(norm_path)
        if not m_sample:
            raise ValueError(f"Could not extract sample id from path:\n{path}")
        sample = int(m_sample.group(1))
        samples.append(sample)

        # --- anglid ---
        m_ang = anglid_pattern.search(norm_path)
        if not m_ang:
            raise ValueError(f"Could not extract anglid from path:\n{path}")
        anglid = int(m_ang.group(1))
        anglids.append(anglid)

        # --- wake_info ---
        # Prefer: folder just above "sampleXXXX"
        wake_info = None
        for i, p in enumerate(parts):
            if p.startswith("sample") and p[6:].isdigit():  # 'sample' + digits
                if i - 1 >= 0:
                    wake_info = parts[i - 1]
                break

        # Fallback: look for a folder starting with '4Mpc_'
        if wake_info is None:
            for p in parts:
                if p.startswith("4Mpc_"):
                    wake_info = p
                    break

        if wake_info is None:
            raise ValueError(f"Could not determine wake_info folder from path:\n{path}")

        wake_infos.append(wake_info)

    return samples, anglids, wake_infos

#%%




def generate_output_paths_void(samples: list[int], anglids: list[int], sliceids: list[int], 
                               wake_infos: list[str], all_data_nowake_void, all_data_wake_void, 
                               range_start: int) -> list[float]:
    """
    Extracts values from preloaded arrays instead of reading files.

    Args:
        samples (list[int]): List of sample numbers.
        anglids (list[int]): List of anglid values (line index).
        sliceids (list[int]): List of slice IDs (position in the line).
        wake_infos (list[str]): List of wake info strings.
        all_data_nowake_void: Preloaded array for no-wake cases.
        all_data_wake_void: Preloaded array for wake cases.
        range_start (int): The first number in the sample range (e.g., 5001 for 5001-5100).

    Returns:
        list[float]: Extracted values from the arrays.
    """

    extracted_values = []

    for sample, anglid, sliceid, wake_info in zip(samples, anglids, sliceids, wake_infos):
        sample_id = sample - range_start  # Generalized index computation

        try:
            if wake_info == '4Mpc_2048c_1024p_zi63_nowakem':
                value = all_data_nowake_void[sample_id, anglid-1, sliceid-1]
            elif wake_info == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m':  # Explicit condition
                value = all_data_wake_void[sample_id, anglid-1, sliceid-1]
            else:
                print(f"Warning: Unrecognized wake_info '{wake_info}' for sample {sample}")
                extracted_values.append(None)
                continue  # Skip to the next iteration

            extracted_values.append(float(value))  # Convert to float
        except IndexError:
            print(f"Warning: Index out of bounds for sample {sample}, anglid {anglid}, sliceid {sliceid}")
            extracted_values.append(None)  # Handle out-of-bounds errors
        except Exception as e:
            print(f"Error accessing data: {e}")
            extracted_values.append(None)  # Handle any other unexpected errors

    return extracted_values


def extract_stat(samples: list[int], anglids: list[int], sliceids: list[int], 
                               wake_infos: list[str], all_data_nowake_void, all_data_wake_void, 
                               rang: list[int]) -> list[float]:
    """
    Extracts values from preloaded arrays instead of reading files.

    Args:
        samples (list[int]): List of sample numbers.
        anglids (list[int]): List of anglid values (line index).
        sliceids (list[int]): List of slice IDs (position in the line).
        wake_infos (list[str]): List of wake info strings.
        all_data_nowake_void: Preloaded array for no-wake cases.
        all_data_wake_void: Preloaded array for wake cases.
        range_start (int): The first number in the sample range (e.g., 5001 for 5001-5100).

    Returns:
        list[float]: Extracted values from the arrays.
    """

    sample_to_idx = {s: i for i, s in enumerate(rang)}

    extracted_values = []

    for sample, anglid, sliceid, wake_info in zip(samples, anglids, sliceids, wake_infos):
        sample_id = sample_to_idx[sample]  # 0-based index into arrays # Generalized index computation

        try:
            if wake_info == '4Mpc_2048c_1024p_zi63_nowakem':
                value = all_data_nowake_void[sample_id, anglid-1, sliceid-1]
            elif wake_info == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m':  # Explicit condition
                value = all_data_wake_void[sample_id, anglid-1, sliceid-1]
            else:
                print(f"Warning: Unrecognized wake_info '{wake_info}' for sample {sample}")
                extracted_values.append(None)
                continue  # Skip to the next iteration

            extracted_values.append(float(value))  # Convert to float
        except IndexError:
            print(f"Warning: Index out of bounds for sample {sample}, anglid {anglid}, sliceid {sliceid}")
            extracted_values.append(None)  # Handle out-of-bounds errors
        except Exception as e:
            print(f"Error accessing data: {e}")
            extracted_values.append(None)  # Handle any other unexpected errors

    return extracted_values

def extract_stat3d(samples: list[int], anglids: list[int], 
                               wake_infos: list[str], all_data_nowake_void, all_data_wake_void, 
                               rang: list[int]) -> list[float]:
    """
    Extracts values from preloaded arrays instead of reading files.

    Args:
        samples (list[int]): List of sample numbers.
        anglids (list[int]): List of anglid values (line index).
        wake_infos (list[str]): List of wake info strings.
        all_data_nowake_void: Preloaded array for no-wake cases.
        all_data_wake_void: Preloaded array for wake cases.
        range_start (int): The first number in the sample range (e.g., 5001 for 5001-5100).

    Returns:
        list[float]: Extracted values from the arrays.
    """

    sample_to_idx = {s: i for i, s in enumerate(rang)}

    extracted_values = []

    for sample, anglid, wake_info in zip(samples, anglids, wake_infos):
        # Skip samples not in rang
        if sample not in sample_to_idx:
            print(f"Skipping sample {sample} (not in rang)")
            continue
        
        sample_id = sample_to_idx[sample]  # 0-based index into arrays # Generalized index computation

        try:
            if wake_info == '4Mpc_2048c_1024p_zi63_nowakem':
                value = all_data_nowake_void[sample_id, anglid-1]
            elif wake_info == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m':  # Explicit condition
                value = all_data_wake_void[sample_id, anglid-1]
            else:
                print(f"Warning: Unrecognized wake_info '{wake_info}' for sample {sample}")
                extracted_values.append(None)
                continue  # Skip to the next iteration

            extracted_values.append(float(value))  # Convert to float
        except IndexError:
            print(f"Warning: Index out of bounds for sample {sample}, anglid {anglid}")
            extracted_values.append(None)  # Handle out-of-bounds errors
        except Exception as e:
            print(f"Error accessing data: {e}")
            extracted_values.append(None)  # Handle any other unexpected errors

    return extracted_values








def extract_stat_diff(samples: list[int], anglids: list[int], sliceids: list[int], 
                               wake_infos: list[str], all_data_nowake_void, all_data_wake_void, 
                               rang: list[int]) -> list[float]:
    """
    Extracts values from preloaded arrays instead of reading files.

    Args:
        samples (list[int]): List of sample numbers.
        anglids (list[int]): List of anglid values (line index).
        sliceids (list[int]): List of slice IDs (position in the line).
        wake_infos (list[str]): List of wake info strings.
        all_data_nowake_void: Preloaded array for no-wake cases.
        all_data_wake_void: Preloaded array for wake cases.
        range_start (int): The first number in the sample range (e.g., 5001 for 5001-5100).

    Returns:
        list[float]: Extracted values from the arrays.
    """
    
    sample_to_idx = {s: i for i, s in enumerate(rang)}


    extracted_values = []

    for sample, anglid, sliceid, wake_info in zip(samples, anglids, sliceids, wake_infos):
        # sample_id = sample - range_start  # Generalized index computation
        sample_id = sample_to_idx[sample]  # 0-based index into arrays # Generalized index computation

        try:
            if wake_info == '4Mpc_2048c_1024p_zi63_nowakem':
                value = all_data_nowake_void[sample_id, anglid-1, sliceid-1] - all_data_wake_void[sample_id, anglid-1, sliceid-1]
            elif wake_info == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m':  # Explicit condition
                value = all_data_wake_void[sample_id, anglid-1, sliceid-1] - all_data_nowake_void[sample_id, anglid-1, sliceid-1]
            else:
                print(f"Warning: Unrecognized wake_info '{wake_info}' for sample {sample}")
                extracted_values.append(None)
                continue  # Skip to the next iteration

            extracted_values.append(float(value))  # Convert to float
        except IndexError:
            print(f"Warning: Index out of bounds for sample {sample}, anglid {anglid}, sliceid {sliceid}")
            extracted_values.append(None)  # Handle out-of-bounds errors
        except Exception as e:
            print(f"Error accessing data: {e}")
            extracted_values.append(None)  # Handle any other unexpected errors

    return extracted_values

def extract_stat_diff3d(samples: list[int], anglids: list[int],
                               wake_infos: list[str], all_data_nowake_void, all_data_wake_void, 
                               rang: list[int]) -> list[float]:
    """
    Extracts values from preloaded arrays instead of reading files.

    Args:
        samples (list[int]): List of sample numbers.
        anglids (list[int]): List of anglid values (line index).
        wake_infos (list[str]): List of wake info strings.
        all_data_nowake_void: Preloaded array for no-wake cases.
        all_data_wake_void: Preloaded array for wake cases.
        range_start (int): The first number in the sample range (e.g., 5001 for 5001-5100).

    Returns:
        list[float]: Extracted values from the arrays.
    """
    
    sample_to_idx = {s: i for i, s in enumerate(rang)}


    extracted_values = []

    for sample, anglid, wake_info in zip(samples, anglids, wake_infos):
        # sample_id = sample - range_start  # Generalized index computation
        
        # Skip samples not in rang
        if sample not in sample_to_idx:
            print(f"Skipping sample {sample} (not in rang)")
            continue
        
        sample_id = sample_to_idx[sample]  # 0-based index into arrays # Generalized index computation

        try:
            if wake_info == '4Mpc_2048c_1024p_zi63_nowakem':
                value = all_data_nowake_void[sample_id, anglid-1] - all_data_wake_void[sample_id, anglid-1]
            elif wake_info == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m':  # Explicit condition
                value = all_data_wake_void[sample_id, anglid-1] - all_data_nowake_void[sample_id, anglid-1]
            else:
                print(f"Warning: Unrecognized wake_info '{wake_info}' for sample {sample}")
                extracted_values.append(None)
                continue  # Skip to the next iteration

            extracted_values.append(float(value))  # Convert to float
        except IndexError:
            print(f"Warning: Index out of bounds for sample {sample}, anglid {anglid}")
            extracted_values.append(None)  # Handle out-of-bounds errors
        except Exception as e:
            print(f"Error accessing data: {e}")
            extracted_values.append(None)  # Handle any other unexpected errors

    return extracted_values


#%%


# def select_extreme_files(data_signal_diff: list[float], wake_infos: list[str], 
#                          files_list: list[str], top_percentage: float = 10) -> tuple[list[str], list[int]]:
#     """
#     Selects the highest X% of positive values where wake is present and the lowest NS negative values where wake is absent.

#     Args:
#         data_signal_diff (list[float]): List of signal differences.
#         wake_infos (list[str]): List indicating wake presence.
#         files_list (list[str]): List of file names.
#         top_percentage (float): Percentage of top positive values to select from the wake group.

#     Returns:
#         tuple[list[str], list[int]]: Selected file names and their positions in files_list.
#     """
    
#     import numpy as np


#     # Convert to numpy arrays for easier indexing
#     data_signal_diff = np.array(data_signal_diff)
#     wake_infos = np.array(wake_infos)
#     files_list = np.array(files_list)

#     # Masks for wake and no-wake conditions
#     wake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
#     nowake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_nowakem'

#     # Filter positive signal differences for wake cases
#     wake_indices = np.where(wake_mask & (data_signal_diff > 0))[0]
#     if len(wake_indices) > 0:
#         wake_values = data_signal_diff[wake_indices]
#         num_select = max(1, int(len(wake_values) * (top_percentage / 100)))  # Select top X%
#         top_wake_indices = wake_indices[np.argsort(wake_values)[-num_select:]]  # Top X% highest
#     else:
#         top_wake_indices = np.array([])

#     # Number of selected items (NS)
#     NS = len(top_wake_indices)

#     # Filter negative signal differences for no-wake cases
#     nowake_indices = np.where(nowake_mask & (data_signal_diff < 0))[0]
#     if len(nowake_indices) > 0:
#         nowake_values = data_signal_diff[nowake_indices]
#         num_nowake = min(NS, len(nowake_values))  # Select at most NS values
#         bottom_nowake_indices = nowake_indices[np.argsort(nowake_values)[:num_nowake]]  # NS lowest
#     else:
#         bottom_nowake_indices = np.array([])

#     # Adjust NS if there are not enough nowake values
#     if len(bottom_nowake_indices) < NS:
#         top_wake_indices = top_wake_indices[:len(bottom_nowake_indices)]  # Reduce wake selection

#     # Combine selected indices
#     selected_indices = np.concatenate([top_wake_indices, bottom_nowake_indices])
    
#     # Get file names, positions and diff values
#     selected_files = files_list[selected_indices].tolist()
#     selected_positions = selected_indices.tolist()
#     selected_signal_diff = data_signal_diff[selected_indices].tolist()
    
#     return selected_files, selected_positions, selected_signal_diff



# # #How to Use:

# # data_signal_diff = [0.5, -0.3, 1.2, -0.8, 0.9, -0.1, 1.5, -1.1, 2.0, -0.5]
# # wake_infos = [
# #     "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m", "4Mpc_2048c_1024p_zi63_nowakem",
# #     "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m", "4Mpc_2048c_1024p_zi63_nowakem",
# #     "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m", "4Mpc_2048c_1024p_zi63_nowakem",
# #     "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m", "4Mpc_2048c_1024p_zi63_nowakem",
# #     "4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m", "4Mpc_2048c_1024p_zi63_nowakem"
# # ]
# # files_list = [f"file_{i}.txt" for i in range(10)]

# # selected_files, selected_positions, selected_signal_diff = select_extreme_files(data_signal_diff, wake_infos, files_list, top_percentage=40)

# # print("Selected Files:", selected_files)
# # print("Selected Positions:", selected_positions)


# import numpy as np

# def select_extreme_files(
#     data_signal_diff: list[float], wake_infos: list[str], 
#     files_list: list[str], data_void: list[float],
#     wake_top_percentage: float = 10, void_percentage: float = 0
# ) -> tuple[list[str], list[int], list[float]]:
#     """
#     Selects:
#       - The top X% highest positive values of data_signal_diff where wake is present.
#       - The lowest negative values of data_signal_diff where wake is absent, after void filtering.
#       - Ensures presence of:
#         - ONLY data_void == 0 entries if void_percentage == 0.
#         - data_void == 0 and the X% smallest values where data_void > 0 if void_percentage > 0.
#       - Adjusts number of wake selections if there are fewer no-wake values.

#     Args:
#         data_signal_diff (list[float]): List of signal differences.
#         wake_infos (list[str]): List indicating wake presence.
#         files_list (list[str]): List of file names.
#         data_void (list[float]): List of void data.
#         wake_top_percentage (float): Percentage of top positive wake values to select.
#         void_percentage (float): Percentage of smallest data_void > 0 to include.

#     Returns:
#         tuple[list[str], list[int], list[float]]: Selected file names, their positions, and corresponding signal differences.
#     """

#     # Convert to numpy arrays for easier indexing
#     data_signal_diff = np.array(data_signal_diff)
#     wake_infos = np.array(wake_infos)
#     files_list = np.array(files_list)
#     data_void = np.array(data_void)

#     # Masks for wake and no-wake conditions
#     wake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
#     nowake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_nowakem'

#     ###### Wake selection ######
#     # Select top X% highest positive wake values
#     wake_indices = np.where(wake_mask & (data_signal_diff > 0))[0]
#     if len(wake_indices) > 0:
#         wake_values = data_signal_diff[wake_indices]
#         num_wake_select = max(1, int(len(wake_values) * (wake_top_percentage / 100)))  # Select top X%
#         top_wake_indices = wake_indices[np.argsort(wake_values)[-num_wake_select:]]  # Top X% highest
#     else:
#         top_wake_indices = np.array([], dtype=int)

#     # Apply void condition for wake
#     void_zero_indices_wake = top_wake_indices[data_void[top_wake_indices] == 0.0]

#     if void_percentage > 0:
#         # Include X% smallest values of data_void > 0 in the wake selection
#         wake_positive_void_indices = top_wake_indices[data_void[top_wake_indices] > 0]
#         num_void_wake_select = max(1, int(len(wake_positive_void_indices) * (void_percentage / 100)))
#         smallest_void_wake_indices = wake_positive_void_indices[np.argsort(data_void[wake_positive_void_indices])[:num_void_wake_select]]
#         top_wake_indices = np.concatenate([void_zero_indices_wake, smallest_void_wake_indices])
#     else:
#         top_wake_indices = void_zero_indices_wake
        


#     num_selected_pairs = len(top_wake_indices)

#     ###### No-wake selection ######
#     # First, apply the void condition before selecting lowest negative values
#     nowake_indices = np.where(nowake_mask & (data_signal_diff < 0))[0]

#     # Filter only data_void == 0
#     void_zero_indices_nowake = nowake_indices[data_void[nowake_indices] == 0.0]

#     if void_percentage > 0:
#         # Include X% smallest values of data_void > 0
#         nowake_positive_void_indices = nowake_indices[data_void[nowake_indices] > 0]
#         num_void_nowake_select = max(1, int(len(nowake_positive_void_indices) * (void_percentage / 100)))
#         smallest_void_nowake_indices = nowake_positive_void_indices[np.argsort(data_void[nowake_positive_void_indices])[:num_void_nowake_select]]
#         nowake_indices_filtered = np.concatenate([void_zero_indices_nowake, smallest_void_nowake_indices])
#     else:
#         nowake_indices_filtered = void_zero_indices_nowake

#     # Now select lowest negative no-wake values, trying to match num_selected_pairs
#     if len(nowake_indices_filtered) > 0:
#         nowake_values = data_signal_diff[nowake_indices_filtered]
#         num_nowake_select = min(num_selected_pairs, len(nowake_values))  # Try to match num_selected_pairs
#         bottom_nowake_indices = nowake_indices_filtered[np.argsort(nowake_values)[:num_nowake_select]]
#     else:
#         bottom_nowake_indices = np.array([], dtype=int)

#     # Adjust number of wake selections if there are fewer no-wake values
#     if len(bottom_nowake_indices) < num_selected_pairs:
#         top_wake_indices = top_wake_indices[:len(bottom_nowake_indices)]  # Reduce wake selection

#     ###### Final selection ######
#     selected_indices = np.concatenate([top_wake_indices, bottom_nowake_indices])
#     selected_files = files_list[selected_indices].tolist()
#     selected_positions = selected_indices.tolist()
#     selected_signal_diff = data_signal_diff[selected_indices].tolist()

#     return selected_files, selected_positions, selected_signal_diff


#%%

import numpy as np

def select_extreme_files(
    data_signal_diff: list[float], wake_infos: list[str], 
    files_list: list[str], data_void: list[float],
    wake_top_percentage: float = 10, void_percentage: float = 0
) -> tuple[list[str], list[int], list[float]]:
    """
    Selects:
      - The top X% highest positive values of data_signal_diff where wake is present.
      - The lowest negative values of data_signal_diff where wake is absent, after void filtering.
      - Ensures presence of:
        - ONLY data_void == 0 entries if void_percentage == 0.
        - data_void == 0 and the X% smallest values where data_void > 0 if void_percentage > 0.
      - Adjusts number of wake selections if there are fewer no-wake values.

    Args:
        data_signal_diff (list[float]): List of signal differences.
        wake_infos (list[str]): List indicating wake presence.
        files_list (list[str]): List of file names.
        data_void (list[float]): List of void data.
        wake_top_percentage (float): Percentage of top positive wake values to select.
        void_percentage (float): Percentage of smallest data_void > 0 to include.

    Returns:
        tuple[list[str], list[int], list[float]]: Selected file names, their positions, and corresponding signal differences.
    """

    # Convert to numpy arrays for easier indexing
    data_signal_diff = np.array(data_signal_diff)
    wake_infos = np.array(wake_infos)
    files_list = np.array(files_list)
    data_void = np.array(data_void)

    # Masks for wake and no-wake conditions
    wake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
    nowake_mask = wake_infos == '4Mpc_2048c_1024p_zi63_nowakem'

    ###### Wake selection ######
    # Select top X% highest positive wake values
    
    # Selection logic:
    # - Wake: keep the strongest positive signal-difference cases
    # - No-wake: keep the most negative cases
    # - Optionally allow a controlled fraction of small nonzero void contamination
    # - Finally balance wake and no-wake counts
    
    
    wake_indices = np.where(wake_mask & (data_signal_diff > 0))[0]

    if len(wake_indices) > 1:
        wake_values = data_signal_diff[wake_indices]
        num_wake_select = max(2, int(len(wake_values) * (wake_top_percentage / 100)))  # Select top X%
        top_wake_indices = wake_indices[np.argsort(wake_values)[-num_wake_select:]]  # Top X% highest
    else:
        wake_indices = np.where(wake_mask)[0]
        wake_values = data_signal_diff[wake_indices]
        top_wake_indices = wake_indices[np.argsort(wake_values)[-2:]]
        


    # Apply void condition for wake
    # void_zero_indices_wake = top_wake_indices[data_void[top_wake_indices] == 0.0]
    
    # print(f"Number of selected elements:={len(top_wake_indices)}")

    
    void_zero_indices_wake = top_wake_indices[np.isclose(data_void[top_wake_indices], 0.0, atol=1e-12)]
    
    # print(f"Number of selected elements:={len(void_zero_indices_wake)}")


    if void_percentage > 0:
        # Include X% smallest values of data_void > 0 in the wake selection
        wake_positive_void_indices = top_wake_indices[data_void[top_wake_indices] > 0]
        num_void_wake_select = max(1, int(len(wake_positive_void_indices) * (void_percentage / 100)))
        smallest_void_wake_indices = wake_positive_void_indices[np.argsort(data_void[wake_positive_void_indices])[:num_void_wake_select]]
        if len(void_zero_indices_wake) + len(smallest_void_wake_indices) > 1:
            top_wake_indices = np.concatenate([void_zero_indices_wake, smallest_void_wake_indices])
    else:
        if len(void_zero_indices_wake) > 1:
            top_wake_indices = void_zero_indices_wake
        
    
        

    # print(f"Number of selected elements:={len(top_wake_indices)}")



    num_selected_pairs = len(top_wake_indices)
    


    ###### No-wake selection ######
    # First, apply the void condition before selecting lowest negative values
    nowake_indices = np.where(nowake_mask & (data_signal_diff < 0))[0]
    

    
    if len(nowake_indices) < 2:
        nowake_indices = np.where(nowake_mask)[0]
        nowake_values = data_signal_diff[nowake_indices]
        nowake_indices = nowake_indices[np.argsort(nowake_values)[-2:]]
        


    
    # print(f"Number of selected elements:={len(nowake_indices)}")


    # Filter only data_void == 0
    void_zero_indices_nowake = nowake_indices[data_void[nowake_indices] == 0.0]
    


    if void_percentage > 0:
        # Include X% smallest values of data_void > 0
        nowake_positive_void_indices = nowake_indices[data_void[nowake_indices] > 0]
        num_void_nowake_select = max(1, int(len(nowake_positive_void_indices) * (void_percentage / 100)))
        smallest_void_nowake_indices = nowake_positive_void_indices[np.argsort(data_void[nowake_positive_void_indices])[:num_void_nowake_select]]
        if len(void_zero_indices_nowake) + len(smallest_void_nowake_indices) > 1:
            nowake_indices_filtered = np.concatenate([void_zero_indices_nowake, smallest_void_nowake_indices])
        else:
            nowake_indices_filtered = nowake_indices
    else:
        if len(void_zero_indices_nowake) > 1:
            nowake_indices_filtered = void_zero_indices_nowake
        else:
            nowake_indices_filtered = nowake_indices
            


    # Now select lowest negative no-wake values, trying to match num_selected_pairs
    if len(nowake_indices_filtered) > 0:
        nowake_values = data_signal_diff[nowake_indices_filtered]
        num_nowake_select = min(num_selected_pairs, len(nowake_values))  # Try to match num_selected_pairs
        bottom_nowake_indices = nowake_indices_filtered[np.argsort(nowake_values)[:num_nowake_select]]
    else:
        bottom_nowake_indices = np.array([], dtype=int)

    # At this point, selected_indices contains a balanced set of
    # difficult/representative wake and no-wake slice candidates.

    # Adjust number of wake selections if there are fewer no-wake values
    if len(bottom_nowake_indices) < num_selected_pairs:
        top_wake_indices = top_wake_indices[:len(bottom_nowake_indices)]  # Reduce wake selection
        
    print(f"Number of selected elements:={len(bottom_nowake_indices)},{len(top_wake_indices)}")


    



    ###### Final selection ######
    selected_indices = np.concatenate([top_wake_indices, bottom_nowake_indices])
    selected_files = files_list[selected_indices].tolist()
    selected_positions = selected_indices.tolist()
    selected_signal_diff = data_signal_diff[selected_indices].tolist()

    return selected_files, selected_positions, selected_signal_diff


# Example
# # Example data
# data_signal_diff = np.array([0.8, 0.5, 0.9, -0.3, -0.2, -0.4, 0.1, -0.1, 0.7, -0.5])
# wake_infos = np.array([
#     '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m',  # Wake
#     '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m',  # Wake
#     '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m',  # Wake
#     '4Mpc_2048c_1024p_zi63_nowakem',             # No-wake
#     '4Mpc_2048c_1024p_zi63_nowakem',             # No-wake
#     '4Mpc_2048c_1024p_zi63_nowakem',             # No-wake
#     '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m',  # Wake
#     '4Mpc_2048c_1024p_zi63_nowakem',             # No-wake
#     '4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m',  # Wake
#     '4Mpc_2048c_1024p_zi63_nowakem'              # No-wake
# ])
# files_list = np.array(["file1", "file2", "file3", "file4", "file5", "file6", "file7", "file8", "file9", "file10"])
# data_void = np.array([0.0, 0.0, 0.2, 0.1, 0.0, 0.3, 0.0, 0.0, 0.4, 0.5])  # Some void contamination

# # Select top 50% wake values, keep void_percentage = 50%
# selected_files, selected_positions, selected_signal_diff = select_extreme_files(
#     data_signal_diff.tolist(), wake_infos.tolist(), files_list.tolist(), data_void.tolist(),
#     wake_top_percentage=50, void_percentage=50
# )

# print("Selected Files:", selected_files)
# print("Selected Positions:", selected_positions)
# print("Selected Signal Differences:", selected_signal_diff)


#%%


def find_extreme_files(selected_signal_diff, selected_files):
    
    import numpy as np


    # Convert to numpy arrays for easier handling
    selected_signal_diff = np.array(selected_signal_diff)
    selected_files = np.array(selected_files)

    # Find the smallest positive value and corresponding file
    positive_mask = selected_signal_diff > 0
    if np.any(positive_mask):
        smallest_positive_index = np.argmin(selected_signal_diff[positive_mask])
        smallest_positive_value = selected_signal_diff[positive_mask][smallest_positive_index]
        smallest_positive_file = selected_files[positive_mask][smallest_positive_index]
    else:
        smallest_positive_value, smallest_positive_file = None, None

    # Find the highest negative value and corresponding file
    negative_mask = selected_signal_diff < 0
    if np.any(negative_mask):
        highest_negative_index = np.argmax(selected_signal_diff[negative_mask])  # Closest to zero
        highest_negative_value = selected_signal_diff[negative_mask][highest_negative_index]
        highest_negative_file = selected_files[negative_mask][highest_negative_index]
    else:
        highest_negative_value, highest_negative_file = None, None

    # Print results
    print(f"Smallest Positive: {smallest_positive_value}, File: {smallest_positive_file}")
    print(f"Highest Negative: {highest_negative_value}, File: {highest_negative_file}")

    return (smallest_positive_value, smallest_positive_file), (highest_negative_value, highest_negative_file)

# # Example usage:
# selected_signal_diff = [0.5, -0.3, 1.2, -0.8, 0.9, -0.1, 1.5, -1.1, 2.0, -0.5]
# selected_files = [f"file_{i}.txt" for i in range(len(selected_signal_diff))]

# find_extreme_files(selected_signal_diff, selected_files)


#%%


# # def split_unique_samples(samples: list[int], files_list: list[str], validation_fraction: float) -> tuple:
# def split_unique_samples(samples: list[int], validation_fraction: float) -> tuple:
#     """
#     Splits unique samples into validation and train_test sets, then retrieves corresponding file names.

#     Args:
#         samples (list[int]): List of sample IDs.
#         files_list (list[str]): List of file names corresponding to samples.
#         validation_fraction (float): Fraction of unique samples for validation.

#     Returns:
#         tuple: (validation_indices, train_test_indices, files_list_validation, files_list_trainTest)
#     """

#     import numpy as np

#     # Get unique sample values
#     unique_samples = list(set(samples))

#     # Determine number of validation samples (at least 1)
#     num_validation = max(1, int(len(unique_samples) * validation_fraction))

#     # Shuffle unique samples for randomness
#     np.random.shuffle(unique_samples)

#     # Split into validation and train_test sets
#     validation_samples = unique_samples[:num_validation]
#     train_test_samples = unique_samples[num_validation:]

#     print(f"Validation Samples ({len(validation_samples)}): {validation_samples}")
#     print(f"Train_Test Samples ({len(train_test_samples)}): {train_test_samples}")

#     # Get indices corresponding to validation and train_test samples
#     validation_indices = [i for i, s in enumerate(samples) if s in validation_samples]
#     train_test_indices = [i for i, s in enumerate(samples) if s in train_test_samples]

#     # # Retrieve file lists
#     # files_list_validation = [files_list[i] for i in validation_indices]
#     # files_list_trainTest = [files_list[i] for i in train_test_indices]

#     # print(f"Validation Indices: {validation_indices}")
#     # print(f"Train_Test Indices: {train_test_indices}")
    
#     # return validation_indices, train_test_indices, files_list_validation, files_list_trainTest
#     # return files_list_validation, files_list_trainTest
#     return validation_indices, train_test_indices

# def split_unique_samples(samples: list[int], wake_infos_all: list[str], validation_fraction: float) -> tuple:
#     """
#     Splits samples into validation and train_test sets, but only keeps samples
#     that contain both wake and no-wake entries.

#     Args:
#         samples (list[int]): List of sample IDs.
#         wake_infos_all (list[str]): Wake labels for each entry.
#         validation_fraction (float): Fraction of unique samples for validation.

#     Returns:
#         tuple: (validation_indices, train_test_indices)
#     """

#     import numpy as np
#     from collections import defaultdict

#     # Group wake labels by sample
#     sample_wake_map = defaultdict(set)

#     for s, w in zip(samples, wake_infos_all):
#         sample_wake_map[s].add(w)

#     # Keep only samples that contain BOTH wake and no-wake
#     unique_samples = [
#         s for s, wakes in sample_wake_map.items()
#         if len(wakes) >= 2
#     ]

#     if len(unique_samples) == 0:
#         raise ValueError("No samples contain both wake and no-wake entries.")

#     # Determine number of validation samples
#     num_validation = max(1, int(len(unique_samples) * validation_fraction))

#     # Shuffle for randomness
#     np.random.shuffle(unique_samples)

#     # Split
#     validation_samples = unique_samples[:num_validation]
#     train_test_samples = unique_samples[num_validation:]

#     print(f"Validation Samples ({len(validation_samples)}): {validation_samples}")
#     print(f"Train_Test Samples ({len(train_test_samples)}): {train_test_samples}")

#     # Retrieve indices
#     validation_indices = [i for i, s in enumerate(samples) if s in validation_samples]
#     train_test_indices = [i for i, s in enumerate(samples) if s in train_test_samples]

#     return validation_indices, train_test_indices

def split_unique_samples(
    samples: list[int],
    wake_infos_all: list[str],
    validation_fraction: float,
    test_fraction: float,
    seed: int = 42,
) -> tuple:
    """
    Splits unique samples into train, validation, and test sets.

    Only keeps samples that contain both wake and no-wake entries.

    Args:
        samples: List of sample IDs.
        wake_infos_all: Wake labels for each entry.
        validation_fraction: Fraction of unique samples for validation.
        test_fraction: Fraction of unique samples for test.
        seed: Random seed for reproducibility.

    Returns:
        tuple: (train_indices, validation_indices, test_indices)
    """

    import numpy as np
    from collections import defaultdict

    # Group wake labels by sample
    sample_wake_map = defaultdict(set)

    for s, w in zip(samples, wake_infos_all):
        sample_wake_map[s].add(w)

    # Keep only samples that contain BOTH wake and no-wake
    unique_samples = [
        s for s, wakes in sample_wake_map.items()
        if len(wakes) >= 2
    ]

    if len(unique_samples) == 0:
        raise ValueError("No samples contain both wake and no-wake entries.")

    if validation_fraction + test_fraction >= 1.0:
        raise ValueError("validation_fraction + test_fraction must be < 1.0")

    # Reproducible shuffle
    rng = np.random.RandomState(seed)
    rng.shuffle(unique_samples)

    n_total = len(unique_samples)

    num_validation = max(1, int(n_total * validation_fraction))
    num_test = max(1, int(n_total * test_fraction))

    test_samples = unique_samples[:num_test]
    validation_samples = unique_samples[num_test:num_test + num_validation]
    train_samples = unique_samples[num_test + num_validation:]

    if len(train_samples) == 0:
        raise ValueError(
            "No training samples left. Reduce validation_fraction or test_fraction."
        )

    print(f"Train Samples ({len(train_samples)}): {train_samples}")
    print(f"Validation Samples ({len(validation_samples)}): {validation_samples}")
    print(f"Test Samples ({len(test_samples)}): {test_samples}")

    train_indices = [i for i, s in enumerate(samples) if s in train_samples]
    validation_indices = [i for i, s in enumerate(samples) if s in validation_samples]
    test_indices = [i for i, s in enumerate(samples) if s in test_samples]

    return train_indices, validation_indices, test_indices

# # Example usage
# samples = [1, 2, 3, 2, 4, 5, 6, 3, 7, 8, 9, 10, 1, 6, 7, 5]
# files_list = [f"file_{i}" for i in range(len(samples))]
# validation_fraction = 0.3

# files_list_validation, files_list_trainTest = split_unique_samples(samples, files_list, validation_fraction)

# print("Validation Files:", val_files)
# print("Train_Test Files:", train_test_files)

#%%






# #%%


# # General Data


# # rang=range(3001,3010+1)
# rang=range(5001,5100+1)
# n_angle = 96


# wake_spec = ["4Mpc_2048c_1024p_zi63_nowakem/","4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/"]


# # Void Data

# # File path and specifications
# path_data = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs/"
# path_void = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_stat/void/"
# path_WakeSignal = "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_stat_2dc1l1_3dc1l1/"


# slices_void = 33
# slices_signal = 32


# percentage_positiveWakeSig = 40

# validation_fraction = 0.5



# #%%

# all_data_nowake_void, all_data_wake_void = data_void(rang, n_angle, slices_void, wake_spec, path_void)
# all_data_nowake_signal, all_data_wake_signal = data_signal(rang, n_angle, slices_signal, wake_spec, path_WakeSignal)




# #%%

# # # simple example

# # # Example usage
# # paths = [
# #     "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs/4Mpc_2048c_1024p_zi63_nowakem/sample5004/sample5004-anglid_1-2dproj_z3_ts32_sl32.png",
# #     "/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5001/sample5001-anglid_86-2dproj_z3_ts32_sl21.png"
# # ]


# # samples, anglids, tilesizes, sliceids, wake_infos = extract_info(paths)
# # print(samples)   # [5004, 5001]
# # print(anglids)   # [1, 86]
# # print(tilesizes) # [32, 32]
# # print(sliceids)  # [32, 21]
# # print(wake_infos)

# # range_start = rang[0]
# # output_paths_void2 = generate_output_paths_void(samples, anglids, sliceids, wake_infos, all_data_nowake_void, all_data_wake_void, range_start)
# # print(output_paths_void2)

# # output_paths_void = extract_stat(samples, anglids, sliceids, wake_infos, all_data_nowake_void, all_data_wake_void, range_start)
# # print(output_paths_void)

# # output_paths_signal= extract_stat(samples, anglids, sliceids, wake_infos, all_data_nowake_signal, all_data_wake_signal, range_start)
# # print(output_paths_signal)

# #%%


# files_list_all = list_all_files(path_data)

# samples_all, anglids_all, tilesizes_all, sliceids_all, wake_infos_all = extract_info(files_list_all)

# range_start = rang[0]

# data_void = extract_stat(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_void, all_data_wake_void, range_start)

# data_signal = extract_stat(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, range_start)

# data_signal_diff = extract_stat_diff(samples_all, anglids_all, sliceids_all, wake_infos_all, all_data_nowake_signal, all_data_wake_signal, range_start)


# #%%

# # split val and trainTest

# files_list_validation, files_list_trainTest = split_unique_samples(samples_all, files_list_all, validation_fraction)


# samples_val, anglids_val, tilesizes_val, sliceids_val, wake_infos_val = extract_info(files_list_validation)
# range_start = rang[0]
# data_void_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_void, all_data_wake_void, range_start)
# data_signal_val = extract_stat(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, range_start)
# data_signal_diff_val = extract_stat_diff(samples_val, anglids_val, sliceids_val, wake_infos_val, all_data_nowake_signal, all_data_wake_signal, range_start)


# samples_tt, anglids_tt, tilesizes_tt, sliceids_tt, wake_infos_tt = extract_info(files_list_trainTest)
# range_start = rang[0]
# data_void_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_void, all_data_wake_void, range_start)
# data_signal_tt = extract_stat(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, range_start)
# data_signal_diff_tt = extract_stat_diff(samples_tt, anglids_tt, sliceids_tt, wake_infos_tt, all_data_nowake_signal, all_data_wake_signal, range_start)



# #%%

# # test plots

# # #  individual values


# # i = 3

# # print(files_list[i])
# # print(data_void[i])
# # print(data_signal[i])




# # # Create an scatter plot with the values, only for the stored data figures.

# # import matplotlib.pyplot as plt
# # import numpy as np

# # # Example data

# # # labels = wake_infos
# # labels = reversed(wake_infos)



# # # Define colors based on labels
# # # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]
# # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]


# # # Scatter plot
# # plt.figure(figsize=(8,6))
# # plt.scatter(data_void, data_signal, c=colors, alpha=0.1, edgecolors='k')

# # # Labels and title
# # plt.xlabel("Data Void")
# # plt.ylabel("Data Signal")
# # plt.title("Scatter Plot of Data Void vs. Data Signal")

# # # Custom legend
# # import matplotlib.patches as mpatches
# # legend_patches = [mpatches.Patch(color='blue', label='No Wake'),
# #                   mpatches.Patch(color='red', label='With Wake')]
# # plt.legend(handles=legend_patches)

# # plt.show()



# # # Create an scatter plot with the values, only for the stored data figures.

# # import matplotlib.pyplot as plt
# # import numpy as np

# # # Example data

# # # labels = wake_infos
# # # labels = reversed(wake_infos)



# # # Define colors based on labels
# # # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]
# # # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]


# # # Scatter plot
# # plt.figure(figsize=(8,6))
# # plt.scatter(all_data_nowake_void[:,:,0:-1].flatten(), all_data_nowake_signal.flatten(), c='blue', alpha=0.1, edgecolors='k')
# # plt.scatter(all_data_wake_void[:,:,0:-1].flatten(), all_data_wake_signal.flatten(), c='red', alpha=0.1, edgecolors='k')

# # # Labels and title
# # plt.xlabel("Data Void")
# # plt.ylabel("Data Signal")
# # plt.title("Scatter Plot of Data Void vs. Data Signal")

# # ax = plt.gca()
# # # ax.set_xlim([xmin, xmax])
# # ax.set_ylim([0, 75])

# # # Custom legend
# # import matplotlib.patches as mpatches
# # legend_patches = [mpatches.Patch(color='blue', label='No Wake'),
# #                   mpatches.Patch(color='red', label='With Wake')]
# # plt.legend(handles=legend_patches)

# # plt.show()







# # # Create an scatter plot with the values signal differences, only for the stored data figures.

# # import matplotlib.pyplot as plt
# # import numpy as np

# # # Example data

# # # labels = wake_infos
# # # labels = reversed(wake_infos)



# # # Define colors based on labels
# # # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]
# # # colors = ['blue' if label == '4Mpc_2048c_1024p_zi63_nowakem' else 'red' for label in labels]


# # # Scatter plot
# # plt.figure(figsize=(8,6))
# # plt.scatter(all_data_wake_void[:,:,0:-1].flatten(), all_data_wake_signal.flatten()-all_data_nowake_signal.flatten(), c='blue', alpha=0.1, edgecolors='k')
# # # plt.scatter(all_data_wake_void[:,:,0:-1].flatten(), all_data_wake_signal.flatten(), c='red', alpha=0.1, edgecolors='k')

# # # Labels and title
# # plt.xlabel("Data Void")
# # plt.ylabel("Data Signal")
# # plt.title("Scatter Plot of Data Void vs. Data Signal")

# # ax = plt.gca()
# # # ax.set_xlim([xmin, xmax])
# # ax.set_ylim([-50, 50])

# # # Custom legend
# # import matplotlib.patches as mpatches
# # legend_patches = [mpatches.Patch(color='blue', label='No Wake'),
# #                   mpatches.Patch(color='red', label='With Wake')]
# # plt.legend(handles=legend_patches)

# # plt.show()


# #%%


# # data_signal_diff

# # wake_infos

# # files_list



# # selected_files, selected_positions, selected_signal_diff = select_extreme_files(data_signal_diff, wake_infos, files_list, percentage_positiveWakeSig)

# # find_extreme_files(selected_signal_diff, selected_files)


# # Select top 50% wake values, keep void_percentage = 50%
# selected_files, selected_positions, selected_signal_diff = select_extreme_files(
#     data_signal_diff, wake_infos_all, files_list_all, data_void,
#     wake_top_percentage=20, void_percentage=0
# )


# find_extreme_files(selected_signal_diff, selected_files)




# #%%

# # validation

# # Select top 50% wake values, keep void_percentage = 50%
# selected_files_val, selected_positions_val, selected_signal_diff_val = select_extreme_files(
#     data_signal_diff_val, wake_infos_val, files_list_validation, data_void_val,
#     wake_top_percentage=20, void_percentage=0
# )


# find_extreme_files(selected_signal_diff_val, selected_files_val)

# # train and test

# # Select top 50% wake values, keep void_percentage = 50%
# selected_files_tt, selected_positions_tt, selected_signal_diff_tt = select_extreme_files(
#     data_signal_diff_tt, wake_infos_tt, files_list_trainTest, data_void_tt,
#     wake_top_percentage=20, void_percentage=0
# )


# find_extreme_files(selected_signal_diff_tt, selected_files_tt)


