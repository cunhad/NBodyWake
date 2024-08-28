#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:22:44 2024

@author: asus
"""



def read_pids(nfiles,ncells,filepath,redshift):
    
    import numpy as np

    
    i = 0
    filename = filepath+redshift+'PID'+str(i)+'.dat'
    node = int(i)
    number_node_dim = nfiles**(1./3)    
    print(filename,node,number_node_dim)
    
    data = np.fromfile(filename, dtype=np.integer, offset=6*(8))
    
    return data

def read_pids_xvs(nfiles,ncells,filepath,redshift):
    
    import numpy as np

    
    i = 1
    filename = filepath+redshift+'PID'+str(i)+'.dat'
 
    # print(filename,node,number_node_dim)
    
    data_pid = np.fromfile(filename, dtype=np.integer, offset=6*(8))
    
    
    filename = filepath+redshift+'xv'+str(i)+'.dat'
    node = int(i)
    nc = ncells
    number_node_dim = nfiles**(1./3)   
    k_node = np.floor(node/number_node_dim**2)
    res = np.floor(node % number_node_dim**2)
    j_node = np.floor(res/number_node_dim);
    i_node=res % number_node_dim
    
    data_xv = np.fromfile(filename, dtype=np.float32 , offset=4*(12)).reshape((-1,6))
    
    data_xv[:,0] = data_xv[:,0] + (nc/number_node_dim)*i_node
    data_xv[:,1] = data_xv[:,1] + (nc/number_node_dim)*j_node
    data_xv[:,2] = data_xv[:,2] + (nc/number_node_dim)*k_node
    
    
    
    return data_pid, data_xv[:,0:3]

def extract_xv_wake(Nmesh,nfiles,ncells,filepath,redshift,npart):
    
    import numpy as np


    grid_spacing = ncells/Nmesh

    # # find the ranges in which the particles in the wake are
    # find the cubes in which the particles of the wake are

    # nfiles = 8
    # npart = 4

    cubes_per_dim = int(round(nfiles ** (1/3)))
    cubes_per_plane = cubes_per_dim * cubes_per_dim
    cubes_below_wake = list(range(int(nfiles / 2 - cubes_per_plane), int(nfiles / 2)))
    cubes_above_wake = list(range(int(nfiles / 2), int(nfiles / 2 + cubes_per_plane)))

    part_wake_id_start = []
    part_wake_id_end = []

    npart_dim_cube = int(npart / cubes_per_dim)
    nun_part_slabcube = int(npart_dim_cube * npart_dim_cube)
    nun_part_cube  = int(npart_dim_cube * npart_dim_cube * npart_dim_cube)

    for cub in cubes_below_wake:
        part_wake_id_start.append(((cub + 1) * nun_part_cube) - nun_part_slabcube + 1)
        part_wake_id_end.append((cub + 1) * nun_part_cube)

    for cub in cubes_above_wake:
        part_wake_id_start.append((cub * nun_part_cube)  + 1)
        part_wake_id_end.append((cub * nun_part_cube) + nun_part_slabcube)

    # 

    # pid_wake_tot = []
    # xv_wake_tot = []

    # Initialize empty arrays to store combined results
    # pid_wake_tot = np.empty((0, 1))  # Start with an empty array of shape (0, 1)
    pid_wake_tot = np.array([])     # Start with an empty vector (1D array)
    xv_wake_tot = np.empty((0, 3))  # Start with an empty array of shape (0, 3)


    # nod = 1
    # for i in range(nod,nod+1):
    for i in range(0,nfiles):    
        
        
        filename = filepath+redshift+'PID'+str(i)+'.dat'        
        # filename = filepath+'PID'+str(i)+'.ic'        
        data_pid = np.fromfile(filename, dtype=np.integer, offset=12*4)
        # data_pid = np.fromfile(filename, dtype=np.integer, offset=4)    #for ic

       
        filename = filepath+redshift+'xv'+str(i)+'.dat'
        # filename = filepath+'xv'+str(i)+'.ic'
        node = int(i)
        nc = ncells
        number_node_dim = nfiles**(1./3)   
        k_node = np.floor(node/number_node_dim**2)
        res = np.floor(node % number_node_dim**2)
        j_node = np.floor(res/number_node_dim);
        i_node=res % number_node_dim
        
        data_xv = np.fromfile(filename, dtype=np.float32 , offset=4*(12)).reshape((-1,6))
        # data_xv = np.fromfile(filename, dtype=np.float32 , offset=4).reshape((-1,6)) # for ic
        
        data_xv[:,0] = (data_xv[:,0] + (nc/number_node_dim)*i_node) / grid_spacing
        data_xv[:,1] = (data_xv[:,1] + (nc/number_node_dim)*j_node) / grid_spacing
        data_xv[:,2] = (data_xv[:,2] + (nc/number_node_dim)*k_node) / grid_spacing
        
        # range_wake = range(np * np * (np - 1) / 2, np * np * (np + 1) / 2)
        
        # List to store index of elements that are within any of the ranges
        within_range_elements = []
        
        # Iterate over each element in lst2
        for index, value in enumerate(data_pid):
            # Check if the element is within any range
            for i in range(len(part_wake_id_start)):
                if part_wake_id_start[i] <= value <= part_wake_id_end[i]:
                    within_range_elements.append(index)
                    break  # Break once we find a range, no need to check further
        
        # indexes_in_range = [index for index, value in enumerate(data_pid) if (npart * npart * (npart - 1) / 2) <= value <= (npart * npart * (npart + 1) / 2)]

        pid_wake = data_pid[within_range_elements]
        xv_wake = data_xv[within_range_elements,0:3]
        
        # pid_wake_tot.extend(pid_wake)
        # xv_wake_tot.extend(* xv_wake)
        
        # Combine with the previous arrays
        pid_wake_tot = np.concatenate([pid_wake_tot, pid_wake])  # Combines n x 1 arrays into ? x 1
        xv_wake_tot = np.vstack([xv_wake_tot, xv_wake])  # Combines n x 3 arrays into ? x 3


        # xv_wake_tot = [combined_list_3xN_row + new_row for combined_list_3xN_row, new_row in zip(xv_wake_tot, xv_wake)]



    
    return xv_wake_tot




def plot_2d_proj_wake(mesh,xv_wake,save=None):
    
    import matplotlib
    # matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    matplotlib.use('Qt5Agg')    # to show figures on desktop
    
    
    from matplotlib import pyplot as plt
    import numpy as np
    
    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    plt.imshow(np.log10(mesh.preview(axes=[1,2])))
    plt.title('2d projection (cell units)')
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    
    # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
    overlay = np.zeros((mesh.shape[2], mesh.shape[1], 4))
    
    alpha=0.3
    # Set the pixels at the positions in pos to red (1, 0, 0)
    for position in xv_wake:
        x, y = int(round(position[2]) % mesh.shape[2]), int(round(position[1]) %  mesh.shape[1])
        overlay[y, x] = [1, 0, 0,alpha]  # Red color
        
    # Overlay the red pixels with 50% transparency
    plt.imshow(overlay)
    
    # # Overlay the positions on the plot
    # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency

    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
        
        return np.log10(mesh.preview(axes=[1,2]))

def obtain_wake_2dgrid_points(mesh,xv_wake):
    
    import numpy as np

    # Apply rounding, modulus, and conversion to int
    grid_points = np.unique(np.vstack([
        np.round(xv_wake[:, 1] % mesh.shape[1]).astype(int),
        np.round(xv_wake[:, 2] % mesh.shape[2]).astype(int)
    ]).T, axis=0)
    
    return grid_points

def obtain_wake_grid_points(mesh,x_wake):
    
    import numpy as np

    # Apply rounding, modulus, and conversion to int
    grid_points_wake = np.unique(np.vstack([
        (np.round(x_wake[:, 0]) % mesh.shape[0]).astype(int),
        (np.round(x_wake[:, 1]) % mesh.shape[1]).astype(int),
        (np.round(x_wake[:, 2]) % mesh.shape[2]).astype(int)
    ]).T, axis=0)
    
    return grid_points_wake



def fraction_of_colapsed(mesh, grid_points_wake, dc_thresh = 1):
    
    import numpy as np

    vals = mesh[grid_points_wake[:,0],grid_points_wake[:,1],grid_points_wake[:,2]]

    tot_grid_wake = len(vals)
    tot_wake_collaps = np.sum(vals >= dc_thresh)
    frac_collaps = tot_wake_collaps / tot_grid_wake
    
    return frac_collaps

def plot_2d_proj_wake_colInfo(mesh,pos_wake,save=None):
    
    import matplotlib
    # matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    matplotlib.use('Qt5Agg')    # to show figures on desktop
    
    grid_points_wake = obtain_wake_grid_points(mesh,pos_wake)
    
    frac_collaps = fraction_of_colapsed(mesh, grid_points_wake)
    
    
    from matplotlib import pyplot as plt
    import numpy as np
    
    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    plt.imshow(np.log10(mesh.preview(axes=[1,2])))
    plt.title(f'2d projection, nc = {frac_collaps:.2f}')
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    
    # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
    overlay = np.zeros((mesh.shape[2], mesh.shape[1], 4))
    
    alpha=0.3
    # Set the pixels at the positions in pos to red (1, 0, 0)
    for position in pos_wake:
        x, y = int(round(position[2]) % mesh.shape[2]), int(round(position[1]) %  mesh.shape[1])
        overlay[y, x] = [1, 0, 0,alpha]  # Red color
        
    # Overlay the red pixels with 50% transparency
    plt.imshow(overlay)
    
    # # Overlay the positions on the plot
    # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency

    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
        
        return np.log10(mesh.preview(axes=[1,2]))


# def den_cont_wake(mesh,grid_points_wake):
    
#     vals = 
    
#     return






