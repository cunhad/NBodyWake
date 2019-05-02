#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 12:09:26 2018

@author: Disrael Cunha
"""

#(run as) python powerSpectrum /home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/ /home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_nowakem/sample0001/ps5/ 32 16
import matplotlib
matplotlib.use('agg')
import sys
import os
import glob
import matplotlib.pyplot as plt
from scipy import interpolate as intp

print(os.path.join(os.path.dirname(__file__), "../../../../Analysis/preprocessing"))

sys.path.append(os.path.join(os.path.dirname(__file__), "../../../Analysis/preprocessing"))

import preprocessing as pps
#file_in=glob.glob('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/snapshot'+'*')
#file_in=glob.glob(sys.argv[1]+sys.argv[2]+'*')

#for x in pps.sorted_nicely(file_in):
#    print(x)

import numpy as np
from pygadgetreader import *

def ReadPos(dir,file):
    pos=readsnap(dir+file,'pos','dm')
    
    return pos

def part2dens3d(part_pos, box_l, bin_x=128):
    """
    Calculate 3D matter density using numpy histograms
    :param part_pos: particle positions in the shape (N, D), where N is particle number and D is dimension
    :param box_l: box length in comoving Mpc/h
    :param bin_x: desired bins per axis for the histogram
    :return: density field
    """
    hist, _edges = np.histogramdd(np.vstack((part_pos[:, 0], part_pos[:, 1], part_pos[:, 2])).T,
                                  bins=bin_x, range=[[0, box_l], [0, box_l], [0, box_l]])
    del _edges
    return hist

def dens2overdens(density, mean_density=None):
    """
    Calculate the overdensity corresponding to a density field
    :param density: input density field
    :param mean_density: if defined normalisation is calculated according to (density - mean(density)) / mean_density
    :return: overdensity field
    """
    assert np.ndim(density) == 3, 'density is not 3D'
    if mean_density:
        delta = (density - np.mean(density)) / mean_density
    else:
        delta = density / np.mean(density) - 1.
    return delta

def power_spectrum(field_x, box_l, bin_k):
    """
        Measures the mass power spectrum of a 3D input field for a given number of bins in Fourier space.
        :param field_x: 3D input field to compute the power spectrum of (typically the overdensity field), dimensionless
        :param box_l: box length of image/cube/box or whatever, units of Mpc or Mpc/h
        :param bin_k: number of bins in Fourier space
        :return: power_k, k: 1D mass power spectrum of field_x, same units as [box_l]**3 and corresponding k values
        """
    assert np.ndim(field_x) == 3, 'field_x is not 3D'
    box_pix = np.size(field_x, axis=0)  # pixel number per axis
    box_dim = np.ndim(field_x)  # dimension

    # This first 'paragraph' is to create masks of indices corresponding to one Fourier bin each
    _freq = np.fft.fftfreq(n=box_pix, d=box_l/box_pix) * 2*np.pi
    _rfreq = np.fft.rfftfreq(n=box_pix, d=box_l/box_pix) * 2*np.pi
    _kx, _ky, _kz = np.meshgrid(_freq, _freq, _rfreq, indexing='ij')
    del _freq, _rfreq
    _k_abs = np.sqrt(_kx**2. + _ky**2. + _kz**2.)
    del _kx, _ky, _kz
    # The following complicated line is actually only creating a 1D array spanning k-space logarithmically from minimum _k_abs to maximum.
    # To start slightly below the minimum and finish slightly above the maximum I use ceil and floor.
    # To ceil and floor not to the next integer but to the next 15th digit, I multiply by 1e15 before flooring and divide afterwards.
    # Since the ceiled/floored value is actually the exponent used for the logspace, going to the next integer would be way too much.
    _k_log = np.logspace(np.floor(np.log10(np.min(_k_abs[1:]))*1.e15)/1.e15, np.ceil(np.log10(np.max(_k_abs[1:]))*1.e15)/1.e15, bin_k)

    _field_k = np.fft.rfftn(np.fft.fftshift(field_x)) * (box_l/box_pix)**box_dim
    power_k = np.empty(np.size(_k_log)-1)
    for i in xrange(np.size(_k_log)-1):
        mask = (_k_abs >= _k_log[i]) & (_k_abs < _k_log[i+1])
        power_k[i] = np.mean(np.abs(_field_k[mask])**2.) / box_l**box_dim
    del _k_abs, _field_k

    k = (_k_log[1:] + _k_log[:-1]) / 2
    del _k_log

    return power_k, k


def plot_ps_gadget(path_in,file_in,path_out,bin_x_,bin_k):
    
    print(path_in+'gadget_out/')
    pos=ReadPos(path_in+'gadget_out/',file_in)
    z=readheader(path_in+'gadget_out/'+file_in,'redshift')
    L=readheader(path_in+'gadget_out/'+file_in,'boxsize')    
    
#    pos=ReadPos('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/','out_z31p000')
#    z=readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/out_z31p000','redshift')
    
    density=part2dens3d(pos, box_l=L, bin_x=bin_x_)
    delta = dens2overdens(density, np.mean(density))
    P_k = power_spectrum(delta, np.float64(L), bin_k)[0]
    k = power_spectrum(delta,  np.float64(L), bin_k)[1]
            
    plt.plot(k,P_k)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel('P(k) [$Mpc^3$]')
    plt.title('Gadget P(k), z=%2.f' %z)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    if not os.path.exists(os.path.dirname(path_out+'gadget/')):
        os.makedirs(path_out+'gadget/')
    plt.savefig(path_out+'gadget/ps_'+file_in+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    return

def plot_ps_picola(path_in,file_in,path_out,bin_x_,bin_k):
    
    print(path_in+'picola_out/')
    pos=ReadPos(path_in+'picola_out/',file_in)
    z=readheader(path_in+'picola_out/'+file_in,'redshift')
    L=readheader(path_in+'picola_out/'+file_in,'boxsize')    
    
#    pos=ReadPos('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/','out_z31p000')
#    z=readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/out_z31p000','redshift')
    
    density=part2dens3d(pos, box_l=L, bin_x=bin_x_)
    delta = dens2overdens(density, np.mean(density))
    P_k = power_spectrum(delta, np.float64(L), bin_k)[0]
    k = power_spectrum(delta,  np.float64(L), bin_k)[1]
            
    plt.plot(k,P_k)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel('P(k) [$Mpc^3$]')
    plt.title('Picola P(k), z=%2.f' %z)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    if not os.path.exists(os.path.dirname(path_out+'picola/')):
        os.makedirs(path_out+'picola/')
    plt.savefig(path_out+'picola/ps_'+file_in+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    return

def plot_ps_comp(path_in,file_in_gadget,file_in_picola,path_out,bin_x_,bin_k):
    
    print(path_in+'gadget_out/')
    pos1=ReadPos(path_in+'gadget_out/',file_in_gadget)
    z=readheader(path_in+'gadget_out/'+file_in_gadget,'redshift')
    L=readheader(path_in+'gadget_out/'+file_in_gadget,'boxsize')  

        
    print(path_in+'picola_out/')
    pos2=ReadPos(path_in+'picola_out/',file_in_picola)
#    z=readheader(path_in+'picola_out/'+file_in_picola,'redshift')
#    L=readheader(path_in+'picola_out/'+file_in_picola,'boxsize') 
    
    density1=part2dens3d(pos1, box_l=L, bin_x=bin_x_)
    delta1 = dens2overdens(density1, np.mean(density1))
    P_k1 = power_spectrum(delta1, np.float64(L), bin_k)[0]
    k1 = power_spectrum(delta1,  np.float64(L), bin_k)[1]
    
    density2=part2dens3d(pos2,box_l=L, bin_x=bin_x_)
    delta2 = dens2overdens(density2, np.mean(density2))
    P_k2 = power_spectrum(delta2, np.float64(L), bin_k)[0]
    k2 = power_spectrum(delta2,  np.float64(L), bin_k)[1]
    
    
    
    plt.plot(k1,P_k1,label="Gadget")
    plt.plot(k2,P_k2,label="MG-Picola")
    plt.legend(loc=1)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel('P(k) [$Mpc^3$]')
    plt.title('Gadget P(k), z=%2.f' %z)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    if not os.path.exists(os.path.dirname(path_out+'pic_gad_comp/')):
        os.makedirs(path_out+'pic_gad_comp/')
    plt.savefig(path_out+'pic_gad_comp/ps_'+file_in_picola+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    plt.plot(k1,P_k1/P_k2)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel(r'$P_{gadget}(k)/P_{picola}(k)$ [$Mpc^3$]')
    plt.title('Frac err of Gadget and Picola P(k), z=%2.f' %z)
    plt.xscale('log')
#    plt.yscale('log')
    #plt.show()
    if not os.path.exists(os.path.dirname(path_out+'pic_gad_comp_ferr/')):
        os.makedirs(path_out+'pic_gad_comp_ferr/')
    plt.savefig(path_out+'pic_gad_comp_ferr/ps_ferr_'+file_in_picola+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    return

def plot_ps_CAMB(redshift):
    
    print('../../../CAMB/transfer_functions/camb_matterpower_z'+str(redshift)+'00.dat')
    P_k=np.loadtxt('../../../CAMB/transfer_functions/camb_matterpower_z'+str(redshift)+'00.dat')
    plt.plot(P_k[:,0],P_k[:,1])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    return P_k

def plot_ps_gadget_CAMB(path_in,file_in_gadget,path_out,bin_x_,bin_k):
    
    print(path_in+'gadget_out/')
    pos1=ReadPos(path_in+'gadget_out/',file_in_gadget)
    z=readheader(path_in+'gadget_out/'+file_in_gadget,'redshift')
    L=readheader(path_in+'gadget_out/'+file_in_gadget,'boxsize')
    
    density1=part2dens3d(pos1, box_l=L, bin_x=bin_x_)
    delta1 = dens2overdens(density1, np.mean(density1))
    P_k1 = power_spectrum(delta1, np.float64(L), bin_k)[0]
    k1 = power_spectrum(delta1,  np.float64(L), bin_k)[1]
    
#    redshift=round((readheader(sys.argv[1]+'gadget_out/'+file_in_gadget,'redshift')))
    k2=np.loadtxt('../../../CAMB/transfer_functions/camb_matterpower_z'+str(round(z))+'00.dat')[:,0]    
    P_k2=np.loadtxt('../../../CAMB/transfer_functions/camb_matterpower_z'+str(round(z))+'00.dat')[:,1]
    
    
    plt.plot(k1,P_k1,label="Gadget")
    plt.plot(k2,P_k2,label="CAMB")
    plt.xlim(xmax = min(k1), xmin = max(k1))
    plt.legend(loc=1)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel('P(k) [$Mpc^3$]')
    plt.title('Gadget and CAMB P(k), z=%2.f' %z)
    plt.xscale('log')
    plt.yscale('log')
#    plt.show()
    if not os.path.exists(os.path.dirname(path_out+'gad_CAMB_comp/')):
        os.makedirs(path_out+'gad_CAMB_comp/')
    plt.savefig(path_out+'gad_CAMB_comp/ps_'+file_in_gadget+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    f2 = intp.interp1d(k2, P_k2, kind='nearest')
    P_k2_=f2(k1)    
    plt.plot(k1,P_k1/P_k2_)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel(r'$P_{gadget}(k)/P_{CAMB}(k)$ [$Mpc^3$]')
    plt.title('Frac of Gadget and CAMB P(k), z=%2.f' %z)
    plt.xscale('log')
#    plt.yscale('log')
#    plt.show()
    if not os.path.exists(os.path.dirname(path_out+'gad_CAMB_comp_ferr/')):
        os.makedirs(path_out+'gad_CAMB_comp_ferr/')
    plt.savefig(path_out+'gad_CAMB_comp_ferr/ps_ferr_'+file_in_gadget+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
        
    return

def plot_ps_picola_CAMB(path_in,file_in_picola,file_in_gadget,path_out,bin_x_,bin_k):
    
        
    print(path_in+'picola_out/')
    pos1=ReadPos(path_in+'picola_out/',file_in_picola)
    z=readheader(path_in+'gadget_out/'+file_in_gadget,'redshift')
    L=readheader(path_in+'picola_out/'+file_in_picola,'boxsize') 
        
    density1=part2dens3d(pos1,box_l=L, bin_x=bin_x_)
    delta1 = dens2overdens(density1, np.mean(density1))
    P_k1 = power_spectrum(delta1, np.float64(L), bin_k)[0]
    k1 = power_spectrum(delta1,  np.float64(L), bin_k)[1]
    
    
    
    k2=np.loadtxt('../../../CAMB/transfer_functions/camb_matterpower_z'+str(round(z))+'00.dat')[:,0]    
    P_k2=np.loadtxt('../../../CAMB/transfer_functions/camb_matterpower_z'+str(round(z))+'00.dat')[:,1]
    
    
    plt.plot(k1,P_k1,label="MG_picola")
    plt.plot(k2,P_k2,label="CAMB")
    plt.xlim(xmax = min(k1), xmin = max(k1))
    plt.legend(loc=1)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel('P(k) [$Mpc^3$]')
    plt.title('MG_picola and CAMB P(k), z=%2.f' %z)
#    plt.xscale('log')
    plt.yscale('log')
#    plt.show()
    if not os.path.exists(os.path.dirname(path_out+'pic_CAMB_comp/')):
        os.makedirs(path_out+'pic_CAMB_comp/')
    plt.savefig(path_out+'pic_CAMB_comp/ps_'+file_in_picola+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
    
    f2 = intp.interp1d(k2, P_k2, kind='nearest')
    P_k2_=f2(k1)    
    plt.plot(k1,P_k1/P_k2_)
    plt.xlabel('k [$Mpc^{-1}$]')
    plt.ylabel(r'$P_{picola}(k)/P_{CAMB}(k)$ [$Mpc^3$]')
    plt.title('Frac of Picola and CAMB P(k), z=%2.f' %z)
    plt.xscale('log')
#    plt.yscale('log')
#    plt.show()
    if not os.path.exists(os.path.dirname(path_out+'pic_CAMB_comp_ferr/')):
        os.makedirs(path_out+'pic_CAMB_comp_ferr/')
    plt.savefig(path_out+'pic_CAMB_comp_ferr/ps_ferr_'+file_in_picola+'.png', bbox_inches = "tight",dpi=300)
    plt.close()
#    
    
    
    
    return


file_in_picola,just_files_picola=pps.list_files_picola(sys.argv[1])    
for x in (just_files_picola):
    print(x)
#    plot_ps_picola('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/',x,'/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_nowakem/sample0001/ps5/',32,16)
    plot_ps_picola(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))

file_in_gadget,just_files_gadget=pps.list_files_gadget(sys.argv[1])
for x in (just_files_gadget):
    print(x)    
#    plot_ps_gadget('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/',x,'/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_nowakem/sample0001/ps5/',32,16)    
    plot_ps_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
    
#print(len(just_files_gadget),len(just_files_picola))    
for x, y in zip(just_files_gadget, just_files_picola):
    print(x, y)
    print(readheader(sys.argv[1]+'gadget_out/'+x,'redshift'),readheader(sys.argv[1]+'picola_out/'+y,'redshift'))
    plot_ps_comp(sys.argv[1],x,y,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
    

#for x in (just_files_gadget[1:-1]):
for x in (just_files_gadget):
    print(x)
    plot_ps_gadget_CAMB(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
   
for x,y in zip(just_files_picola,just_files_gadget:
    print(x,y)
    plot_ps_picola_CAMB(sys.argv[1],x,y,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
        
    
    