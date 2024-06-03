#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:44:01 2023

@author: asus
"""

def plot_2d_proj(mesh,save=None):
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)

    from matplotlib import pyplot as plt
    import numpy as np
    
    # plt.figure()    
    # plt.imshow(mesh.preview(axes=[0,2]))
    plt.imshow(np.log10(mesh.preview(axes=[0,2])))
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
        
        return np.log10(mesh.preview(axes=[0,2]))
    
def plot_2d_proj_eachSlice(mesh,slice_list,dept,save=None):
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut window forwards (non iteractive, like slurm)
    # matplotlib.use('gplot')
    # matplotlib.use('Qt5Agg')    #show figures on desktop

    from matplotlib import pyplot as plt
    import numpy as np
    
    size0 = mesh.value.shape[0]
    size1 = mesh.value.shape[1]
    size2 = mesh.value.shape[2]
    
    array_3d = mesh.value.reshape(size0, int(size1/dept), dept, size2).sum(axis=2)
    
    if save != None:
        splited = save[0].split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    for i in slice_list:
        values = array_3d[:,i,:]
        # plt.figure()  
        plt.imshow(np.log10(values))
        plt.title('2d projection (cell units)')
        # plt.show()
        if save != None:            
            plt.savefig(save[i], bbox_inches = "tight",dpi=300)
            plt.close()
        
        
    return array_3d[:,slice_list,:]