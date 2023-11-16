#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 11:04:34 2023

@author: asus
"""








def histogram_3d(mesh,save=None):
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
    plt.title('3d histogram')
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