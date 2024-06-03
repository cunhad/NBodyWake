#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 11:04:34 2023

@author: asus
"""








def histogram_3d(mesh,save=None):
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)
    
    import numpy
    from matplotlib import pyplot as plt

    # histogram of 1+delta in log-spaced bins

    # one_plus_delta = mesh.paint(mode='real')
    # one_plus_delta = mesh
    values = mesh.value.flatten()
    
    
    
    num_zero = numpy.count_nonzero(values == 0.)
    #remove zeroes
    values = values[values != 0.]
    # values[~numpy.isin(values, 0.)]
    # minim = numpy.count_nonzero(values == 0.)
    # minim = values == 0.
    # minim = values

    # bins = numpy.logspace(-7, numpy.log10(30.), 100)
    minim_log = numpy.log10(values.min())
    maxim_log = numpy.log10(values.max())
    num_per_log10 = 10
    numb_bins = int((maxim_log-minim_log)*num_per_log10)
    bins = numpy.logspace(minim_log, maxim_log, numb_bins)    
    # print(one_plus_delta.value.min())
    hist3d = plt.hist(values, bins=bins)
    # hist3d = plt.hist(values)

    # format the axes
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$1+\delta$")
    plt.ylabel(r"$N_\mathrm{cells}$")
    # plt.xlim(1e-2, 20)
    plt.title("3d histogram, num of 0: "  + str(num_zero))
    plt.xticks(numpy.power(10,numpy.arange(numpy.floor(minim_log), numpy.ceil(maxim_log)+1, 1)))

    
    
    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()
    
    return hist3d


def histogram_3d_eachSlice(mesh,slice_list,dept,save=None):
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)
     
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
        values = values.flatten()
        num_zero = np.count_nonzero(values == 0.)
        values = values[values != 0.]
        minim_log = np.log10(values.min())
        maxim_log = np.log10(values.max())
        num_per_log10 = 10
        numb_bins = int((maxim_log-minim_log)*num_per_log10)
        bins = np.logspace(minim_log, maxim_log, numb_bins)    
        hist2d = plt.hist(values, bins=bins)
        # format the axes
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r"$1+\delta$")
        plt.ylabel(r"$N_\mathrm{cells}$")
        # plt.xlim(1e-2, 20)
        plt.title("3d histogram, num of 0: "  + str(num_zero))
        plt.xticks(np.power(10,np.arange(np.floor(minim_log), np.ceil(maxim_log)+1, 1)))
        # plt.show()
        if save != None:            
            plt.savefig(save[i], bbox_inches = "tight",dpi=300)
            plt.close()
    
    return
