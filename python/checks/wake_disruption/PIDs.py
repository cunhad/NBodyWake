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






def extract_pid_wake(Nmesh, nfiles, ncells, filepath, redshift, npart):
    
    # import Read_slices
    import numpy as np

    
    Nmesh_x, Nmesh_y, Nmesh_z = Nmesh
    # grid_spacing_x = ncells / Nmesh_x
    # grid_spacing_y = ncells / Nmesh_y
    # grid_spacing_z = ncells / Nmesh_z
    # wake_grid_z_start = np.floor(-1+(Nmesh_z)/2)
    # wake_grid_z_end   = np.ceil(+1+(Nmesh_z)/2)
    wake_grid_z_start = np.floor(-3+(Nmesh_z)/2)
    wake_grid_z_end   = np.ceil(-1+(Nmesh_z)/2)
    
    
    
    pid_wake_tot = []
    x_wake_tot = []
    
    number_node_dim = nfiles ** (1/3)
    node = np.arange(nfiles)
    
    k_node = np.floor(node / number_node_dim**2).astype(int)
    res = node % number_node_dim**2
    j_node = np.floor(res / number_node_dim).astype(int)
    i_node = res % number_node_dim

    for i in range(nfiles):    
        # print(i)
        
        data_pid = np.fromfile(filepath + redshift + f'PID{i}.dat', dtype=np.int64, offset=12*4)
        data_xv = np.fromfile(filepath + redshift + f'xv{i}.dat', dtype=np.float32, offset=4*12).reshape((-1, 6))

        data_xv[:, 0] = (data_xv[:, 0] + (ncells / number_node_dim) * i_node[i]) 
        data_xv[:, 1] = (data_xv[:, 1] + (ncells / number_node_dim) * j_node[i]) 
        data_xv[:, 2] = (data_xv[:, 2] + (ncells / number_node_dim) * k_node[i]) 
        
        # Create a mask for elements that fall within the specified ranges in data_pid
        mask = np.zeros(data_xv[:,2].shape, dtype=bool)
        # for start, end in zip(part_wake_id_start, part_wake_id_end):
        mask |= (data_xv[:,2] >= wake_grid_z_start) & (data_xv[:,2] <= wake_grid_z_end)
        
        # Apply the mask to data_pid and data_xv
        pid_wake = data_pid[mask]
        x_wake = data_xv[mask, :3]  # Only apply the mask to the rows of data_xv

        # Append the results to the lists
        pid_wake_tot.append(pid_wake)
        x_wake_tot.append(x_wake)
        
    pid_wake_tot = np.concatenate(pid_wake_tot)
    x_wake_tot = np.vstack(x_wake_tot)
    
    return pid_wake_tot, x_wake_tot


# pos wake particles in the CUBEP3M grid size
def extract_pos_wake(Nmesh, nfiles, ncells, filepath, redshift, npart, pid_wake):
    
    import numpy as np

    # pid_wake_tot = []
    x_wake_tot = []

    number_node_dim = nfiles ** (1/3)
    node = np.arange(nfiles)
    
    k_node = np.floor(node / number_node_dim**2).astype(int)
    res = node % number_node_dim**2
    j_node = np.floor(res / number_node_dim).astype(int)
    i_node = res % number_node_dim

    for i in range(nfiles):    
        # print(i)
        
        data_pid = np.fromfile(filepath + redshift + f'PID{i}.dat', dtype=np.int64, offset=12*4)
        data_xv = np.fromfile(filepath + redshift + f'xv{i}.dat', dtype=np.float32, offset=4*12).reshape((-1, 6))

        data_xv[:, 0] = (data_xv[:, 0] + (ncells / number_node_dim) * i_node[i]) 
        data_xv[:, 1] = (data_xv[:, 1] + (ncells / number_node_dim) * j_node[i])
        data_xv[:, 2] = (data_xv[:, 2] + (ncells / number_node_dim) * k_node[i]) 
        
        # Create a mask that only picks data_pid values that are in pid_wake
        mask = np.isin(data_pid, pid_wake)
        
        # Apply the mask to data_pid and data_xv
        # pid_wake = data_pid[mask]
        x_wake = data_xv[mask, :3]  # Only apply the mask to the rows of data_xv

        # Append the results to the lists
        # pid_wake_tot.append(pid_wake)
        x_wake_tot.append(x_wake)

    # pid_wake_tot = np.concatenate(pid_wake_tot)
    x_wake_tot = np.vstack(x_wake_tot)

    return x_wake_tot





def plot_2d_proj_wake(mesh,grid_points_wake,save=None):
    
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
    for position in grid_points_wake:
        x, y = int(round(position[2]) % mesh.shape[2]), int(round(position[1]) %  mesh.shape[1])
        overlay[x, y] = [1, 0, 0, alpha]  # Red color
        
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
        
        # return np.log10(mesh.preview(axes=[1,2]))
        return

# def obtain_wake_2dgrid_points(Nmesh,xv_wake):
    
#     import numpy as np

#     # Apply rounding, modulus, and conversion to int
#     grid_points = np.unique(np.vstack([
#         np.round(xv_wake[:, 1] % Nmesh[1]).astype(int),
#         np.round(xv_wake[:, 2] % Nmesh[2]).astype(int)
#     ]).T, axis=0)
    
#     return grid_points



# def obtain_wake_grid_points(Nmesh,x_wake):
    
#     import numpy as np

#     # Apply rounding, modulus, and conversion to int
#     grid_points_wake = np.unique(np.vstack([
#         (np.round(x_wake[:, 0]) % Nmesh[0]).astype(int),
#         (np.round(x_wake[:, 1]) % Nmesh[1]).astype(int),
#         (np.round(x_wake[:, 2]) % Nmesh[2]).astype(int)
#     ]).T, axis=0)
    
#     return grid_points_wake


def obtain_wake_grid_points_chunks(Nmesh, pos, chunk_size=10000):
    
    """
    Return unique mesh cells touched by particles using:
      - CIC in X,Y
      - single-slice binning in Z
      - periodic boundaries in X,Y
    """
    import numpy as np
    
    
    Nx, Ny, Nz = map(int, Nmesh)
    unique_cells = set()
    
    for i in range(0, len(pos), chunk_size):
        chunk = pos[i:i+chunk_size]
    
        # CIC shift (as in your MATLAB)
        x = chunk[:, 0] - 0.5
        y = chunk[:, 1] - 0.5
    
        # Z slice (1-based → 0-based)
        z_slice = np.ceil(chunk[:, 2]).astype(int) - 1
        valid = (z_slice >= 0) & (z_slice < Nz)
    
        if not np.any(valid):
            continue
    
        x = x[valid]
        y = y[valid]
        z = z_slice[valid]
    
        # Base and neighbor cells
        i1 = np.floor(x).astype(int) + 1
        j1 = np.floor(y).astype(int) + 1
        i2 = i1 + 1
        j2 = j1 + 1
    
        # Periodic wrapping
        I1 = np.mod(i1, Nx)
        I2 = np.mod(i2, Nx)
        J1 = np.mod(j1, Ny)
        J2 = np.mod(j2, Ny)
    
        # Collect the 4 CIC cells per particle
        cells = np.stack([
            np.stack([I1, J1, z], axis=1),
            np.stack([I2, J1, z], axis=1),
            np.stack([I1, J2, z], axis=1),
            np.stack([I2, J2, z], axis=1),
        ], axis=1).reshape(-1, 3)
    
        unique_cells.update(map(tuple, cells))

    return np.array(list(unique_cells), dtype=int)

    #previous code#
    
    # # Initialize an empty set to store unique grid points
    # unique_grid_points = set()
    
    # # Process x_wake in chunks
    # for i in range(0, len(x_wake), chunk_size):
    #     chunk = x_wake[i:i+chunk_size]
        
    #     # Compute the grid points for the chunk
    #     grid_points_chunk = np.vstack([
    #         (np.round(chunk[:, 0]) % Nmesh[0]).astype(int),
    #         (np.round(chunk[:, 1]) % Nmesh[1]).astype(int),
    #         (np.round(chunk[:, 2]) % Nmesh[2]).astype(int)
    #     ]).T
        
    #     # Add the unique grid points to the set
    #     unique_grid_points.update(map(tuple, grid_points_chunk))
    
    # # Convert the set of unique grid points back to a numpy array
    # grid_points_wake = np.array(list(unique_grid_points))
    
    # return grid_points_wake
    
    
    
    
    

# def obtain_wake_grid_points_eachslice(Nmesh,x_wake):
    
#     import numpy as np

#     # Initialize an empty list of arrays for each possible value in range(0, mesh.shape[0])
#     grid_points_wake_list = [np.empty((0, 2), dtype=int) for _ in range(Nmesh[0])]


#     # Apply rounding, modulus, and conversion to int
#     grid_points_wake = np.unique(np.vstack([
#         (np.round(x_wake[:, 0]) % Nmesh[0]).astype(int),
#         (np.round(x_wake[:, 1]) % Nmesh[1]).astype(int),
#         (np.round(x_wake[:, 2]) % Nmesh[2]).astype(int)
#     ]).T, axis=0)

#     # Sort the grid points into the corresponding lists based on the first coordinate
#     for point in grid_points_wake:
#         first_coordinate = point[0]
#         grid_points_wake_list[first_coordinate] = np.vstack([grid_points_wake_list[first_coordinate], point[1:3]])



#     return grid_points_wake_list


def obtain_wake_grid_points_eachslice_chunks(Nmesh, x_wake, chunk_size=10000):
    
    import numpy as np

    # Initialize an empty list of sets for each possible value in range(0, Nmesh[0])
    grid_points_wake_list = [set() for _ in range(Nmesh[2])]

    # Process x_wake in chunks
    for i in range(0, len(x_wake), chunk_size):
        chunk = x_wake[i:i+chunk_size]
        
        # Compute the grid points for the chunk
        grid_points_chunk = np.vstack([
            (np.round(chunk[:, 0]) % Nmesh[0]).astype(int),
            (np.round(chunk[:, 1]) % Nmesh[1]).astype(int),
            (np.round(chunk[:, 2]) % Nmesh[2]).astype(int)
        ]).T
        
        # Sort the grid points into the corresponding sets based on the first coordinate
        for point in grid_points_chunk:
            first_coordinate = point[2]
            grid_points_wake_list[first_coordinate].add(tuple(point[0:2]))

    # Convert sets back to numpy arrays for final output, ensuring sorted order for consistency
    grid_points_wake_list = [np.array(sorted(s), dtype=int) for s in grid_points_wake_list]
    
    # # Convert sets back to numpy arrays for final output
    # grid_points_wake_list = [np.array(list(s), dtype=int) for s in grid_points_wake_list]
    
    return grid_points_wake_list


def fraction_of_colapsed(mesh, grid_points_wake, dc_thresh = 0.5):
    
    import numpy as np

    vals = mesh[grid_points_wake[:,0],grid_points_wake[:,1],grid_points_wake[:,2]]

    tot_grid_wake = len(vals)
    tot_wake_collaps = np.sum(vals >= dc_thresh)
    frac_collaps = tot_wake_collaps / tot_grid_wake
    
    return frac_collaps

def density_contrast_inside_wake(mesh, grid_points_wake):
    
    import numpy as np

    
    vals = mesh[grid_points_wake[:,0],grid_points_wake[:,1],grid_points_wake[:,2]]
    
    # # # Clip values in vals to ensure they are <= 1
    # vals = np.clip(vals, None, 1)
    
    dcw = sum(vals)/len(vals)
    
    return dcw

def plot_2d_proj_wake_colInfo3d(mesh,grid_points_wake,save=None):
    
    import matplotlib
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
    else:
        matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    
    # Nmesh = [mesh.shape[0],mesh.shape[1],mesh.shape[2]]
    # grid_points_wake = obtain_wake_grid_points(Nmesh,pos_wake)
    
    frac_collaps = fraction_of_colapsed(mesh, grid_points_wake)
    dcw = density_contrast_inside_wake(mesh, grid_points_wake)
    
    from matplotlib import pyplot as plt
    import numpy as np
    
    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    # plt.imshow(np.log10(+1+mesh.preview(axes=[0,1])))
    plt.imshow(np.log10(+1+np.sum(mesh, axis=2)/mesh.shape[2]))
    
    # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
    plt.title('2d projection, nc3d = {:.2f}, dcw = {:.2f}'.format(frac_collaps, dcw))
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    
    # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
    overlay = np.zeros((mesh.shape[0], mesh.shape[1], 4))
    
    alpha=0.3
    # Set the pixels at the positions in pos to red (1, 0, 0)
    for position in grid_points_wake:
        x, y = int(round(position[0]) % mesh.shape[0]), int(round(position[1]) %  mesh.shape[1])
        overlay[x, y] = [1, 0, 0,alpha]  # Red color
        
    # Overlay the red pixels with 50% transparency
    plt.imshow(overlay)
    
    # # Overlay the positions on the plot
    # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency

    if save is not None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
        
        # return np.log10(mesh.preview(axes=[1,2]))
    
    return


def plot_2d_proj_wakediff_colInfo3d(mesh, mesh_nowake, grid_points_wake,save=None):
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
    else:
        matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    
    # Nmesh = [mesh.shape[0],mesh.shape[1],mesh.shape[2]]
    # grid_points_wake = obtain_wake_grid_points(Nmesh,pos_wake)
    
    dcwd = density_contrast_inside_wake(mesh - mesh_nowake, grid_points_wake)
    
    clipdiff = np.clip(mesh, None, 1) - np.clip(mesh_nowake, None, 1)
    # min_clipdiff = np.min(clipdiff)
    
    dcwdc = density_contrast_inside_wake(clipdiff, grid_points_wake)
    frac_collaps = fraction_of_colapsed(clipdiff, grid_points_wake)

    
    

    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    # plt.imshow(np.log10(+1+mesh.preview(axes=[0,1])))
    # plt.imshow(np.log10(+1+np.sum(mesh, axis=2)/mesh.shape[2]))
    plt.imshow(np.sum(clipdiff, axis=2)/mesh.shape[2])

    
    # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
    plt.title('2d projection, nc3dd = {:.2f}, dcwd = {:.2f}, dcwdc = {:.2f}'.format(frac_collaps,dcwd, dcwdc))
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    
    # # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
    # overlay = np.zeros((mesh.shape[0], mesh.shape[1], 4))
    
    # alpha=0.3
    # # Set the pixels at the positions in pos to red (1, 0, 0)
    # for position in grid_points_wake:
    #     x, y = int(round(position[0]) % mesh.shape[0]), int(round(position[1]) %  mesh.shape[1])
    #     overlay[x, y] = [1, 0, 0,alpha]  # Red color
        
    # # Overlay the red pixels with 50% transparency
    # plt.imshow(overlay)
    
    # # Overlay the positions on the plot
    # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency

    if save is not None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
        
        # return np.log10(mesh.preview(axes=[1,2]))
    
    return


def fraction_of_colapsed_in2d(mesh2d, grid_points_wake2d, dc_thresh = 0.5):
    
    import numpy as np
    
    # print(mesh2d.shape)
    
    # Check if grid_points_wake2d is empty
    if len(grid_points_wake2d) == 0:
        return 0

    if len(grid_points_wake2d) == 3:
        # print("here")
        vals = mesh2d[grid_points_wake2d[0],grid_points_wake2d[1]]
        tot_grid_wake = vals.size
    else:
        # print(grid_points_wake2d)
        vals = mesh2d[grid_points_wake2d[:,0],grid_points_wake2d[:,1]]
        tot_grid_wake = len(vals)

    
    tot_wake_collaps = np.sum(vals >= dc_thresh)
    frac_collaps = tot_wake_collaps / tot_grid_wake
    
    return frac_collaps

def density_contrast_inside_wake_in2d(mesh2d, grid_points_wake2d):
    
    # import numpy as np

    # Check if grid_points_wake2d is empty
    if len(grid_points_wake2d) == 0:
        return 0
    
    vals = mesh2d[grid_points_wake2d[:,0],grid_points_wake2d[:,1]]
    
    # vals = np.clip(vals, None, 1)

    
    dcw = sum(vals)/len(vals)
    
    return dcw


def plot_2d_slice_wake_colInfo2d(mesh2d,grid_points_wake2d,save=None):
    
    import matplotlib
    # matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    matplotlib.use('Qt5Agg')    # to show figures on desktop
    
    Nmesh2d = [mesh2d.shape[0],mesh2d.shape[1]]
    
    frac_collaps = fraction_of_colapsed_in2d(mesh2d, grid_points_wake2d)
    dcw = density_contrast_inside_wake_in2d(mesh2d, grid_points_wake2d)
    
    from matplotlib import pyplot as plt
    import numpy as np
    
    plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    plt.imshow(np.log10(mesh2d))
    # plt.title(f'2d projection, nc2d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
    plt.title('2d projection, nc2d = {:.2f}, dcw = {:.2f}'.format(frac_collaps, dcw))
    # plt.title('2d slice')
    # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
    
    # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
    overlay = np.zeros((mesh2d.shape[1], mesh2d.shape[0], 4))
    
    alpha=0.3
    # Set the pixels at the positions in pos to red (1, 0, 0)
    for position in grid_points_wake2d:
        x, y = int(round(position[1]) % mesh2d.shape[1]), int(round(position[0]) %  mesh2d.shape[0])
        overlay[x, y] = [1, 0, 0, alpha]  # Red color
        
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
        
        # return np.log10(mesh2d)
        return

# def den_cont_wake(mesh,grid_points_wake):
    
#     vals = 
    
#     return




def rotate_pos_wake(Pos, rot_angle, pivot, nc, nup, resol_factor, depth, lenght_factor=1):
    
    import numpy as np
    
    # Ensure pivot is a NumPy array to handle element-wise operations
    pivot = np.array(pivot)
    
    phi, theta, psi  = rot_angle

    Pos = np.mod(Pos, nc)
    Pos = Pos * (nup * resol_factor) / nc

    axis_size = np.array([
        [nup * resol_factor, 0, 0],
        [0, nup * resol_factor, 0],
        [0, 0, nup * resol_factor]
    ])

    # Subtract from Pos, ensuring that pivot is a NumPy array
    Pos -= (nup * resol_factor / 2) + pivot * (nup * resol_factor / nc)

    Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

    R = Rz @ Ry @ Rx
    Pos = Pos @ R.T
    
    axis_size = axis_size @ R.T

    Pos += (1 / (2 * lenght_factor)) * nup * resol_factor

    lim = (1 / lenght_factor) * nup * resol_factor
    Pos_expand = np.empty((0, 3))

    D = np.array([[Dx, Dy, Dz] for Dx in range(-1, 2) for Dy in range(-1, 2) for Dz in range(-1, 2)])
    # D = np.array([[Dx, Dy, Dz] for Dx in range(0, 1) for Dy in range(0, 1) for Dz in range(0, 1)])
    shifts = D @ axis_size

    for shift in shifts:
        Pos_aux = Pos + shift
        mask = np.all((Pos_aux >= 0) & (Pos_aux < lim), axis=1)
        Pos_expand = np.vstack([Pos_expand, Pos_aux[mask]])
        # mask = np.all((Pos_aux >= 0) & (Pos_aux <= lim))
        # Pos_expand = np.vstack([Pos_expand, Pos_aux[mask]])
        # Pos_expand = np.vstack([Pos_expand, Pos_aux])
        
    # Pos_expand[:,0:2] -= 0.5      
    Pos_expand[:,2] = np.floor(Pos_expand[:,2]/(lim/depth))
    
    # Pos[:,2] = Pos[:,2]/(lim/depth)

    # return Pos_expand, shift, lim, Pos
    return Pos_expand

def plot_2d_proj_onlywake_eachSlice(mesh, grid_points_wake, slice_list,save=None):
    
    threshold = 0.5
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
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
        
        if len(slice_list) == 1:
            mesh_2d = mesh.squeeze()
            # grid_points_wake_2d = grid_points_wake.squeeze()
        else:
            mesh_2d = mesh[:,:,i]
        
        # print(mesh_2d.shape)
        
        # grid_points_wake_2d = grid_points_wake[i]
        
        meshaux = np.asarray(mesh_2d)
        
        nonzero = meshaux != 0
        N_nonzero = np.count_nonzero(nonzero)
        
        above = meshaux > threshold
        
        # 1) Fraction of non-zero cells above threshold
        frac_collaps = np.count_nonzero(above & nonzero) / N_nonzero
        
        # 2) Sum of above-threshold values normalized by non-zero count
        dcw = meshaux[above].sum() / N_nonzero
    
        # frac_collaps = fraction_of_colapsed_in2d(mesh_2d, grid_points_wake_2d)
        # dcw = density_contrast_inside_wake_in2d(mesh_2d, grid_points_wake_2d)
        
        
        
        plt.figure()    
        # plt.imshow(mesh.preview(axes=[0,2]))
        plt.imshow(np.minimum(mesh_2d, 1))
        # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
        plt.title('2d projection, nc3d = {:.2f}, dcw = {:.2f}'.format(frac_collaps, dcw))
        # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
        # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
        
        # (no overlay)
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


def plot_2d_proj_wake_colInfo3d_eachSlice(mesh, grid_points_wake, slice_list,save=None):
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
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
        
        if len(slice_list) == 1:
            mesh_2d = mesh.squeeze()
            # grid_points_wake_2d = grid_points_wake.squeeze()
        else:
            mesh_2d = mesh[:,:,i]
        
        # print(mesh_2d.shape)
        
        grid_points_wake_2d = grid_points_wake[i]
    
        frac_collaps = fraction_of_colapsed_in2d(mesh_2d, grid_points_wake_2d)
        dcw = density_contrast_inside_wake_in2d(mesh_2d, grid_points_wake_2d)
        
        
        
        plt.figure()    
        # plt.imshow(mesh.preview(axes=[0,2]))
        plt.imshow(np.log10(+1+mesh_2d))
        # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
        plt.title('2d projection, nc3d = {:.2f}, dcw = {:.2f}'.format(frac_collaps, dcw))
        # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
        # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
        
        # Create an overlay with zeros (same shape as arr, with 3 color channels for RGB), 4th channel is for alpha
        overlay = np.zeros((mesh_2d.shape[0], mesh_2d.shape[1], 4))
        
        alpha=0.3
        # Set the pixels at the positions in pos to red (1, 0, 0)
        for position in grid_points_wake_2d:
            x, y = int(round(position[0]) % mesh_2d.shape[0]), int(round(position[1]) %  mesh_2d.shape[1])
            overlay[x, y] = [1, 0, 0,alpha]  # Red color
            
        # Overlay the red pixels with 50% transparency
        plt.imshow(overlay)
        
        
        
        # # Overlay the positions on the plot
        # plt.scatter(xv_wake[:,2], xv_wake[:,1], color='red', alpha=0.1)  # alpha=0.5 for 50% transparency
    
        if save is not None:            
            plt.savefig(save[i], bbox_inches = "tight",dpi=300)
            plt.close()
        
    return 
        

def plot_2d_proj_wakediff_colInfo3d_eachSlice(mesh,mesh_nowake, grid_points_wake, slice_list,save=None):
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
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
        
        if len(slice_list) == 1:
            mesh_2d = mesh.squeeze()
            mesh_2d_nowake = mesh_nowake.squeeze()
            # grid_points_wake_2d = grid_points_wake.squeeze()
        else:
            mesh_2d = mesh[:,:,i]
            mesh_2d_nowake = mesh_nowake[:,:,i]
        
        # print(mesh_2d.shape)
        
        grid_points_wake_2d = grid_points_wake[i]
    
        # dcw = density_contrast_inside_wake_in2d(mesh_2d, grid_points_wake_2d)
        
        
                # clipdiff = np.clip(mesh_2d, None, 1) - np.clip(mesh_2d_nowake, None, 1)
        # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, None, 1) 
        # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, -1, 1)
        # mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16)
        # min_clipdiff = np.min(clipdiff)
        
        mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16) - (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
        # Get the minimum and maximum of the array
        min_val = np.min(mesh_test)
        max_val = np.max(mesh_test)
        
        # # Shift and rescale the array to be between -1 and 1
        # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
        # Create a copy of the original array to modify
        mesh_rescaled = np.copy(mesh_test)
        
        # Rescale the values >= 0 to the range [0, 1]
        mask_positive = mesh_rescaled >= 0
        mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
        
        # Rescale the values < 0 to the range [-1, 0]
        mask_negative = mesh_rescaled < 0
        mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val


        
        
        frac_collaps = fraction_of_colapsed_in2d(mesh_rescaled, grid_points_wake_2d)
        dcwdc = density_contrast_inside_wake_in2d(mesh_rescaled, grid_points_wake_2d)
        

        mesh_test2 = (2/np.pi)*np.arctan((mesh_2d+1)*16)- (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
        # mesh_test2 = (2/np.pi)*np.arctan((mesh_2d-mesh_2d_nowake)*32)

        
        # Get the minimum and maximum of the array
        min_val = np.min(mesh_test2)
        max_val = np.max(mesh_test2)
        
        # # Shift and rescale the array to be between -1 and 1
        # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
        # Create a copy of the original array to modify
        mesh_rescaled = np.copy(mesh_test2)
        
        # # Rescale the values >= 0 to the range [0, 1]
        # mask_positive = mesh_rescaled >= 0
        # mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
        # # mesh_rescaled[mesh_rescaled > 0] = 1
        
        # # # Set all values < 0 to 0
        # # mesh_rescaled[mesh_rescaled < 0] = 0
        # # Rescale the values < 0 to the range [-1, 0]
        # mask_negative = mesh_rescaled < 0
        # mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val

        
        
        dcwd = density_contrast_inside_wake_in2d(mesh_rescaled, grid_points_wake_2d)

        
        
        plt.figure()    
        # plt.imshow(mesh.preview(axes=[0,2]))
        # plt.imshow(np.log10(+1+mesh_2d))
        img = plt.imshow(mesh_rescaled)
        plt.colorbar(img, orientation='vertical')  # Add a vertical colorbar on the right
        # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
        plt.title('2d projection, nc3dd = {:.2f}, dcwd = {:.2f},  dcwdc = {:.2f}'.format(frac_collaps, dcwd,dcwdc))
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



def plot_2d_proj_curveletFilt_eachSlice(mesh,mesh_nowake, grid_points_wake, slice_list,save=None):
    
    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    
    import sys
    sys.path.append('/home/asus/Programs/PyCurvelab-master/')

    import pyct
    
    # 
    if save is None:
        matplotlib.use('Qt5Agg')    # to show figures on desktop
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
            
    # Initialize the fdct2 object
    nbs = 2   # Number of scales
    nba = 20  # Number of angles at the 2nd coarsest scale
    ac = True # Use curvelets at the coarsest scale
    
    # print(list(slice_list))
    for i ,ls in enumerate(slice_list):
        
        # print(i)
        
        if len(slice_list) == 1:
            mesh_2d = mesh.squeeze()
            mesh_2d_nowake = mesh_nowake.squeeze()
            # grid_points_wake_2d = grid_points_wake.squeeze()
        else:
            mesh_2d = mesh[:,:,i]
            mesh_2d_nowake = mesh_nowake[:,:,i]
        
        # print(mesh_2d.shape)
        
        shape = mesh_2d.shape
        
        grid_points_wake_2d = grid_points_wake[i]
        
        curvelet_transform = pyct.fdct2(shape, nbs, nba, ac, norm=False, vec=True, cpx=False)
        curvelet_coefficients = curvelet_transform.fwd(mesh_2d)
        
        filtered_arr = curvelet_transform.inv(curvelet_coefficients)

    
        dcw = density_contrast_inside_wake_in2d(filtered_arr, grid_points_wake_2d)
        
        
        # # clipdiff = np.clip(mesh_2d, None, 1) - np.clip(mesh_2d_nowake, None, 1)
        # # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, None, 1) 
        # # clipdiff = np.clip(mesh_2d - mesh_2d_nowake, -1, 1)
        # # mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16)
        # # min_clipdiff = np.min(clipdiff)
        
        # mesh_test = (2/np.pi)*np.arctan((mesh_2d+1)*16) - (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
        # # Get the minimum and maximum of the array
        # min_val = np.min(mesh_test)
        # max_val = np.max(mesh_test)
        
        # # # Shift and rescale the array to be between -1 and 1
        # # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
        # # Create a copy of the original array to modify
        # mesh_rescaled = np.copy(mesh_test)
        
        # # Rescale the values >= 0 to the range [0, 1]
        # mask_positive = mesh_rescaled >= 0
        # mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
        
        # # Rescale the values < 0 to the range [-1, 0]
        # mask_negative = mesh_rescaled < 0
        # mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val


        
        
        frac_collaps = fraction_of_colapsed_in2d(filtered_arr, grid_points_wake_2d)
        dcwdc = density_contrast_inside_wake_in2d(filtered_arr, grid_points_wake_2d)
        

        # mesh_test2 = (2/np.pi)*np.arctan((mesh_2d+1)*16)- (2/np.pi)*np.arctan((mesh_2d_nowake+1)*16) 
        
        # # mesh_test2 = (2/np.pi)*np.arctan((mesh_2d-mesh_2d_nowake)*32)

        
        # # Get the minimum and maximum of the array
        # min_val = np.min(mesh_test2)
        # max_val = np.max(mesh_test2)
        
        # # # Shift and rescale the array to be between -1 and 1
        # # mesh_test2 = 2 * (mesh_test2 - min_val) / (max_val - min_val) - 1
        
        # # Create a copy of the original array to modify
        # mesh_rescaled = np.copy(mesh_test2)
        
        # # # Rescale the values >= 0 to the range [0, 1]
        # # mask_positive = mesh_rescaled >= 0
        # # mesh_rescaled[mask_positive] = mesh_rescaled[mask_positive] / max_val
        # # # mesh_rescaled[mesh_rescaled > 0] = 1
        
        # # # # Set all values < 0 to 0
        # # # mesh_rescaled[mesh_rescaled < 0] = 0
        # # # Rescale the values < 0 to the range [-1, 0]
        # # mask_negative = mesh_rescaled < 0
        # # mesh_rescaled[mask_negative] = mesh_rescaled[mask_negative] / -min_val

        
        
        dcwd = density_contrast_inside_wake_in2d(filtered_arr, grid_points_wake_2d)

        
        
        plt.figure()    
        # plt.imshow(mesh.preview(axes=[0,2]))
        # plt.imshow(np.log10(+1+mesh_2d))
        img = plt.imshow(filtered_arr)
        plt.colorbar(img, orientation='vertical')  # Add a vertical colorbar on the right
        # plt.title(f'2d projection, nc3d = {frac_collaps:.2f}, dcw = {dcw:.2f}')
        plt.title('2d projection, nc3dd = {:.2f}, dcwd = {:.2f},  dcwdc = {:.2f}'.format(frac_collaps, dcwd,dcwdc))
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





import numpy as np


def cic_deposit_cell_units(pos, Nmesh, weights=None, dtype=np.float64):
    """
    CIC deposit onto a 3D mesh with periodic boundaries.
    Positions are in *cell units*: x in [0, Nx), y in [0, Ny), z in [0, Nz).

    Parameters
    ----------
    pos : (n,3) float array
        Particle positions in cell units.
    Nmesh : (3,) iterable of int
        (Nx, Ny, Nz)
    weights : None or (n,) float array
        Optional particle weights (mass). If None, weight=1 for each particle.
    dtype : numpy dtype
        Output mesh dtype.

    Returns
    -------
    mesh : (Nx,Ny,Nz) array
        CIC-deposited mesh.
    """
    
    pos = np.asarray(pos, dtype=np.float64)
    Nx, Ny, Nz = map(int, Nmesh)

    n = pos.shape[0]
    if weights is None:
        w = np.ones(n, dtype=np.float64)
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape != (n,):
            raise ValueError("weights must have shape (n,)")

    # Wrap periodic boundaries (handles negatives and pos == Nx, etc.)
    gx = np.mod(pos[:, 0], Nx)
    gy = np.mod(pos[:, 1], Ny)
    gz = np.mod(pos[:, 2], Nz)

    # Base (left/lower/back) cell
    ix0 = np.floor(gx).astype(np.int64)
    iy0 = np.floor(gy).astype(np.int64)
    iz0 = np.floor(gz).astype(np.int64)

    # Fractional part inside the cell
    tx = gx - ix0
    ty = gy - iy0
    tz = gz - iz0

    # Neighbor cell (right/upper/front) with periodic wrap
    ix1 = (ix0 + 1) % Nx
    iy1 = (iy0 + 1) % Ny
    iz1 = (iz0 + 1) % Nz

    # 1D CIC weights
    wx0, wx1 = (1.0 - tx), tx
    wy0, wy1 = (1.0 - ty), ty
    wz0, wz1 = (1.0 - tz), tz

    mesh = np.zeros((Nx, Ny, Nz), dtype=dtype)

    # Deposit into 8 corners
    np.add.at(mesh, (ix0, iy0, iz0), w * wx0 * wy0 * wz0)
    np.add.at(mesh, (ix0, iy0, iz1), w * wx0 * wy0 * wz1)
    np.add.at(mesh, (ix0, iy1, iz0), w * wx0 * wy1 * wz0)
    np.add.at(mesh, (ix0, iy1, iz1), w * wx0 * wy1 * wz1)

    np.add.at(mesh, (ix1, iy0, iz0), w * wx1 * wy0 * wz0)
    np.add.at(mesh, (ix1, iy0, iz1), w * wx1 * wy0 * wz1)
    np.add.at(mesh, (ix1, iy1, iz0), w * wx1 * wy1 * wz0)
    np.add.at(mesh, (ix1, iy1, iz1), w * wx1 * wy1 * wz1)

    return mesh

def cic_xy_with_zslice(pos, Nmesh, weights=None, shift_xy=0.5, dtype=np.float64):
    """
    Match MATLAB scheme:
      x,y: CIC onto 4 neighbors with x = pos_x - 0.5, y = pos_y - 0.5
      z: choose ONE slice via z_slice = ceil(pos_z)  (1..Nz)
      periodic wrap in x,y via mod
      z is NOT periodic; particles outside (0, Nz] are ignored

    Parameters
    ----------
    pos : (n,3) float array
        Positions in cell units. Example: x=2.5 is inside 3rd cell (index 2).
    Nmesh : (3,) iterable of int
        (Nx, Ny, Nz)
    weights : None or (n,) float array
        Optional particle weights. If None, each particle has weight 1.
    shift_xy : float
        The "-0.5" shift used in your MATLAB.
    dtype : numpy dtype

    Returns
    -------
    mesh : (Nx, Ny, Nz) array
    """
    pos = np.asarray(pos, dtype=np.float64)
    Nx, Ny, Nz = map(int, Nmesh)
    n = pos.shape[0]

    if weights is None:
        w = np.ones(n, dtype=np.float64)
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape != (n,):
            raise ValueError("weights must have shape (n,)")

    # MATLAB: x = Pos(:,1)-0.5; y = Pos(:,2)-0.5
    x = pos[:, 0] - shift_xy
    y = pos[:, 1] - shift_xy

    # MATLAB equivalent with lim=Nz and slice=Nz => dz = 1:
    # z_slice = ceil(z)  (1-based)
    z_slice_1based = np.ceil(pos[:, 2]).astype(np.int64)

    # Keep only valid z slices: (z_slice>0 && z_slice<=Nz)
    m = (z_slice_1based > 0) & (z_slice_1based <= Nz)
    if not np.any(m):
        return np.zeros((Nx, Ny, Nz), dtype=dtype)

    x = x[m]; y = y[m]
    z0 = z_slice_1based[m] - 1  # convert to 0-based
    w = w[m]

    # MATLAB: i1=floor(x)+1; j1=floor(y)+1; i2=i1+1; j2=j1+1
    i1 = np.floor(x).astype(np.int64) + 1
    j1 = np.floor(y).astype(np.int64) + 1
    i2 = i1 + 1
    j2 = j1 + 1

    # MATLAB: dx1=i1-x; dy1=j1-y; dx2=1-dx1; dy2=1-dy1
    dx1 = i1 - x
    dy1 = j1 - y
    dx2 = 1.0 - dx1
    dy2 = 1.0 - dy1

    # Periodic wrap in x,y (MATLAB: mod(i,nb)+1)
    I1 = np.mod(i1, Nx)
    I2 = np.mod(i2, Nx)
    J1 = np.mod(j1, Ny)
    J2 = np.mod(j2, Ny)

    mesh = np.zeros((Nx, Ny, Nz), dtype=dtype)

    # Deposit into 4 neighbors in xy, fixed z slice
    np.add.at(mesh, (I1, J1, z0), w * dx1 * dy1)
    np.add.at(mesh, (I2, J1, z0), w * dx2 * dy1)
    np.add.at(mesh, (I1, J2, z0), w * dx1 * dy2)
    np.add.at(mesh, (I2, J2, z0), w * dx2 * dy2)

    return mesh



