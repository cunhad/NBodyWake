#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:44:01 2023

@author: asus
"""

def plot_2d_proj(mesh,save=None):
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)

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
        