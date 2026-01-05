#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 10:35:40 2024

@author: asus
"""

import sys
sys.path.append('/home/asus/Programs/PyCurvelab-master/')

import pyct


import numpy as np


# help(pyct.fdct2)

# # a = pyct.fdct2.fwd(mesh[:,:,0])


arr = np.log10(+1+mesh[:,:,0])

# arr = mesh[:,:,0]


#%%
# Define the shape of the array
shape = arr.shape
m, n = arr.shape

nbs = int(np.floor(np.log2(min(m,n)))-3)
# b = pyct.fdct2.param(a)


# Initialize the fdct2 object
# nbs = 7   # Number of scales
nba = 16  # Number of angles at the 2nd coarsest scale
# ac = True # Use curvelets at the coarsest scale
# curvelet_transform = pyct.fdct2(shape, nbs, nba, ac=False, norm=False, vec=True, cpx=False)
curvelet_transform = pyct.fdct2(shape, nbs+1, nba, ac=False, norm=False, vec=False, cpx=False)

# Apply the forward curvelet transform
curvelet_coefficients = curvelet_transform.fwd(arr)



# filtered_arr = curvelet_transform.inv(curvelet_coefficients)



#%%

scales = len(curvelet_coefficients)
# low_sc = 1
# high_scale = scales
# scale_now = 1


curvelet_coefficients_filt = curvelet_coefficients

for scale in range(0,scales-2):
# for scale in range(scales-2,scales):
    print("scale is ", scale)
    numb_angl = len(curvelet_coefficients[scale])

    for angle in range(0,numb_angl):
        print(angle)
        curvelet_coefficients_filt[scale][angle] = np.zeros_like(curvelet_coefficients[scale][angle])

filtered_arr = curvelet_transform.inv(curvelet_coefficients_filt)






#%%
from matplotlib import pyplot as plt

import matplotlib
matplotlib.use('TkAgg')    # to show figures on desktop

# import os    
# os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
#%%

plt.figure()    
# plt.imshow(mesh.preview(axes=[0,2]))
# plt.imshow(np.log10(+1+mesh_2d))
img = plt.imshow(arr)
plt.colorbar(img, orientation='vertical')  # Add a vertical colorbar on the right
plt.show()
#%%

plt.figure()    
# plt.imshow(mesh.preview(axes=[0,2]))
# plt.imshow(np.log10(+1+mesh_2d))
img = plt.imshow(filtered_arr)
plt.colorbar(img, orientation='vertical')  # Add a vertical colorbar on the right
plt.show()

#%%




# curvelet_coefficients now contains the curvelet transform of 'arr'

#%%

# def plot_2d_proj_curveletFilt_eachSlice(mesh,mesh_nowake, grid_points_wake, slice_list,save=None):
    
#     import matplotlib
#     from matplotlib import pyplot as plt
#     import numpy as np
    
#     # 
#     if save is None:
#         matplotlib.use('Qt5Agg')    # to show figures on desktop
#     else:
#         matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    
#     # Nmesh = [mesh.shape[0],mesh.shape[1],mesh.shape[2]]
#     # grid_points_wake = obtain_wake_grid_points(Nmesh,pos_wake)
    
#     if save != None:
#         splited = save[0].split('/')   
#         folder = "/".join(splited[0:-1])
#         import os
#         if not os.path.exists(folder):
#             os.makedirs(folder)
            
#     # Initialize the fdct2 object
#     nbs = 2   # Number of scales
#     nba = 16  # Number of angles at the 2nd coarsest scale
#     ac = True # Use curvelets at the coarsest scale
    
#     # print(list(slice_list))
#     for i ,ls in enumerate(slice_list):
        
#         # print(i)
        
#         if len(slice_list) == 1:
#             mesh_2d = mesh.squeeze()
#             mesh_2d_nowake = mesh_nowake.squeeze()
#             # grid_points_wake_2d = grid_points_wake.squeeze()
#         else:
#             mesh_2d = mesh[:,:,i]
#             mesh_2d_nowake = mesh_nowake[:,:,i]
        
#         # print(mesh_2d.shape)
        
#         shape = mesh_2d.shape
        
#         grid_points_wake_2d = grid_points_wake[i]
        
#         curvelet_transform = pyct.fdct2(shape, nbs, nba, ac, norm=False, vec=True, cpx=False)
#         curvelet_coefficients = curvelet_transform.fwd(mesh_2d)
        
#         filtered_arr = curvelet_transform.inv(curvelet_coefficients)

    
#         dcw = density_contrast_inside_wake_in2d(filtered_arr, grid_points_wake_2d)
        
        
#         # # clipdiff = np.clip(mesh_2d, None, 1) - np.clip(mesh_2d_nowake, None, 1)
#         # # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, None, 1) 
#         # # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, -1, 1)
#         # # mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16)
#         # # min_clipdiff = np.min(clipdiff)
        
#         # mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16) - (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
#         # # Get the minimum and maximum of the array
#         # min_val = np.min(mesh_test)
#         # max_val = np.max(mesh_test)
        
#         # # # Shift and rescale the array to be between -1 and 1
#         # # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
#         # # Create a copy of the original array to modify
#         # mesh_rescaled = np.copy(mesh_test)
        
#         # # Rescale the values >= 0 to the range [0, 1]
#         # mask_positive = mesh_rescaled >= 0
#         # mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
        
#         # # Rescale the values < 0 to the range [-1, 0]
#         # mask_negative = mesh_rescaled < 0
#         # mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val


        
        
#         frac_collaps = fraction_of_colapsed_in2d(filtered_arr, grid_points_wake_2d)
#         dcwdc = density_contrast_inside_wake_in2d(filtered_arr, grid_points_wake_2d)
        

#         # mesh_test2 = (2/np.pi)*np.arctan((mesh_2d+1)*16)- (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
#         # # mesh_test2 = (2/np.pi)*np.arctan((mesh_2d-mesh_2d_nowake)*32)

        
#         # # Get the minimum and maximum of the array
#         # min_val = np.min(mesh_test2)
#         # max_val = np.max(mesh_test2)
        
#         # # # Shift and rescale the array to be between -1 and 1
#         # # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
#         # # Create a copy of the original array to modify
#         # mesh_rescaled = np.copy(mesh_test2)
        
#         # # # Rescale the values >= 0 to the range [0, 1]
#         # # mask_positive = mesh_rescaled >= 0
#         # # mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
#         # # # mesh_rescaled[mesh_rescaled > 0] = 1
        
#         # # # # Set all values < 0 to 0
#         # # # mesh_rescaled[mesh_rescaled < 0] = 0
#         # # # Rescale the values < 0 to the range [-1, 0]
#         # # mask_negative = mesh_rescaled < 0
#         # # mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val

        
        
#         # dcwd = density_contrast_inside_wake_in2d(mesh_rescaled, grid_points_wake_2d)

        
        
#         plt.figure()    
#         # plt.imshow(mesh.preview(axes=[0,2]))
#         # plt.imshow(np.log10(+1+mesh_2d))
#         img = plt.imshow(mesh_rescaled)
#         plt.colorbar(img, orientation='vertical')  # Add a vertical colorbar on the right
#         # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
#         plt.title('2d projection, nc3dd = {:.2f}, dcwd = {:.2f},  dcwdc = {:.2f}'.format(frac_collaps, dcwd,dcwdc))
#         # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
#         # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
        
#         # # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
#         # overlay = np.zeros((mesh_2d.shape[0], mesh_2d.shape[1], 4))
        
#         # alpha=0.3
#         # # Set the pixels at the positions in pos to red (1, 0, 0)
#         # for position in grid_points_wake_2d:
#         #     x, y = int(round(position[0]) % mesh_2d.shape[0]), int(round(position[1]) %  mesh_2d.shape[1])
#         #     overlay[x, y] = [1, 0, 0,alpha]  # Red color
            
#         # # Overlay the red pixels with 50% transparency
#         # plt.imshow(overlay)
        
        
        
#         # # Overlay the positions on the plot
#         # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency
    
#         if save is not None:            
#             plt.savefig(save[i], bbox_inches = "tight",dpi=300)
#             plt.close()
        
#     return 