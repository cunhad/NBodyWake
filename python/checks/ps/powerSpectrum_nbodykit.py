#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:06:17 2023

@author: Disrael
"""

from nbodykit.lab import *

def powerSpectrum(mesh,save=None):
    # r = FFTPower(mesh, mode='1d', dk=0.05, kmin=0.1)
    r = FFTPower(mesh, mode='1d')
    Pk = r.power

    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)
    
    # #print out the meta-data
    # for k in Pk.attrs:
    #     print("%s = %s" %(k, str(Pk.attrs[k])))
    from matplotlib import pyplot as plt
    
    plt.figure()    
    # plt.loglog(Pk['k'],(pow(Pk['k'],3))*Pk['power'].real)
    # plt.loglog(Pk['k'],Pk['power'].real)
    # or print the shot noise subtracted P(k)
    plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'])
    plt.xlabel(r"$k$ $[h \mathrm{Mpc}^{-1}]$")
    plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^{3}]$")
    # plt.title('Dimensionless Power Spectrum')
    plt.title('Power Spectrum')
    # plt.title('Gadget, Picola, CAMB and EH P(k), z=%2.f' %z)
    
    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)            
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()

    return Pk




def powerSpectrum2d(mesh,nmu,save=None):
    # r = FFTPower(mesh, mode='1d', dk=0.05, kmin=0.1)
    # r = FFTPower(mesh, mode='1d')
    # r = FFTPower(mesh, mode='2d', dk=0.005, kmin=0.01, Nmu=5, los=[0,0,1])
    r = FFTPower(mesh, mode='2d',Nmu=nmu)    
    Pkmu = r.power

    # print(Pkmu)
    # print(Pkmu.coords)
    # # plot each mu bin
    # for i in range(Pkmu.shape[1]):
    #     Pk = Pkmu[:,i] # select the ith mu bin
    #     label = r'$\mu$=%.1f' % (Pkmu.coords['mu'][i])
    #     plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'], label=label)
    
    import matplotlib
    matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)
    
    from matplotlib import pyplot as plt

    plt.figure()
    # plot each mu bin
    for i in range(Pkmu.shape[1]):
        Pk = Pkmu[:,i] # select the ith mu bin
        label = r'$\mu$=%.1f' % (Pkmu.coords['mu'][i])
        plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'], label=label)    
    # format the axes
    plt.legend(loc=0, ncol=2)
    plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    plt.ylabel(r"$P(k, \mu)$ [$h^{-3}\mathrm{Mpc}^3$]")
    plt.xlim(0.01, 0.6)
    
    if save != None:
        splited = save.split('/')   
        folder = "/".join(splited[0:-1])
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)            
        plt.savefig(save, bbox_inches = "tight",dpi=300)
        plt.close()


    return Pkmu