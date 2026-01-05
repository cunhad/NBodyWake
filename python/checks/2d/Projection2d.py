#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:44:01 2023

@author: asus
"""

def plot_2d_proj(mesh,save=None):
    
    import matplotlib
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
    else:
        matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)

    from matplotlib import pyplot as plt
    import numpy as np
    
    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    # plt.imshow(np.log10(1+mesh.preview(axes=[0,1])))
    plt.imshow(np.log10(+1+np.sum(mesh, axis=2)/mesh.shape[2]))
    plt.title('2d projection (cell units)')
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
        

    # return np.log10(1+np.sum(mesh, axis=2))
    return
       
    
# def plot_2d_proj_eachSlice(mesh,slice_list,dept,save=None):
    
#     import matplotlib
#     matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
#     # matplotlib.use('gplot')
#     # matplotlib.use('Qt5Agg')    #show figures on desktop

#     from matplotlib import pyplot as plt
#     import numpy as np
    
#     size0 = mesh.value.shape[0]
#     size1 = mesh.value.shape[1]
#     size2 = mesh.value.shape[2]
    
#     array_3d = mesh.value.reshape(size0, int(size1/dept), dept, size2).sum(axis=2)
    
#     if save != None:
#         splited = save[0].split('/')   
#         folder = "/".join(splited[0:-1])
#         import os
#         if not os.path.exists(folder):
#             os.makedirs(folder)
    
#     for i in slice_list:
#         values = array_3d[:,i,:]
#         # plt.figure()  
#         plt.imshow(np.log10(values))
#         plt.title('2d projection (cell units)')
#         # plt.show()
#         if save != None:            
#             plt.savefig(save[i], bbox_inches = "tight",dpi=300)
#             plt.close()
        
        
#     return array_3d[:,slice_list,:]


def plot_2d_proj_eachSlice(mesh, slice_list,save=None):
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
        # plt.switch_backend('Agg')
    else:        
        matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    
    # Nmesh = [mesh.shape[0],mesh.shape[1],mesh.shape[2]]
    # grid_points_wake = obtain_wake_grid_points(Nmesh,pos_wake)
    
    if save != None:
        splited = save[0].split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    # print(list(slice_list))
    for i ,ls in enumerate(slice_list):
        
        # print(i)
        
        if len(slice_list) == 1 & np.shape(mesh)[2]==1:
            mesh_2d = mesh.squeeze()
            # grid_points_wake_2d = grid_points_wake.squeeze()
        else:
            mesh_2d = mesh[:,:,i]
        
        # print(mesh_2d.shape)
        
        # grid_points_wake_2d = grid_points_wake[i]
    

        
        plt.figure()    
        # plt.imshow(mesh.preview(axes=[0,2]))
        plt.imshow(np.log10(+1+mesh_2d))
        # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
        plt.title('2d projection')
        # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
        # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
        
        # # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
        # overlay = np.zeros((mesh_2d.shape[0], mesh_2d.shape[1], 4))
        
        # alpha=0.3
        # # Set the pixels at the positions in pos to red (1, 0, 0)
        # for position in grid_points_wake_2d:
        #     x, y = int(round(position[0]) % mesh_2d.shape[0]), int(round(position[1]) %  mesh_2d.shape[1])
        #     overlay[x, y] = [1, 0, 0,alpha]  # Red color
            
        # # Overlay the red pixels with 50% transparency
        # plt.imshow(overlay)
        
        
        
        # # Overlay the positions on the plot
        # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency
    
        if save is not None:            
            plt.savefig(save[i], bbox_inches = "tight",dpi=300)
            plt.close()
        
    return 
        