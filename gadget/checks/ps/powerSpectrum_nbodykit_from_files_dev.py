#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 21:01:44 2019

@author: Disrael
"""


import matplotlib
#matplotlib inline
#config InlineBackend.figure_format = 'retina'




#from nbodykit.lab import cosmology

from nbodykit.lab import *
from nbodykit import setup_logging, style

import matplotlib.pyplot as plt
plt.style.use(style.notebook)

#setup_logging()
#
#redshift = 0.55
#cosmo = cosmology.Planck15
#Plin = cosmology.LinearPower(cosmo, redshift, transfer='EisensteinHu')
#b1 = 2.0
#
#cat = LogNormalCatalog(Plin=Plin, nbar=3e-4, BoxSize=1380., Nmesh=128, bias=b1, seed=42)
#
#
## add RSD
#line_of_sight = [0,0,1]
#cat['RSDPosition'] = cat['Position'] + cat['VelocityOffset'] * line_of_sight
#
## convert to a MeshSource, using TSC interpolation on 256^3 mesh
#mesh = cat.to_mesh(window='tsc', Nmesh=256, compensated=True, position='RSDPosition')
#




#import matplotlib
matplotlib.use('agg')
import sys
import os
import glob
#import matplotlib.pyplot as plt
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


from nbodykit.lab import cosmology

cosmo = cosmology.Cosmology()
cosmo = cosmo.clone(h=0.7, nonlinear=True,T0_cmb= 2.7255,Omega0_b=0.0445,Omega0_cdm=0.246,n_s=0.96,sigma8=0.8628)
print("new sigma8 = %.4f" % cosmo.sigma8)










def ReadPos(dir,file):
    pos=readsnap(dir+file,'pos','dm')
    
    return pos

def ReadVel(dir,file):
    vel=readsnap(dir+file,'vel','dm')
    
    return vel

def ReadMass(dir,file):
    mass=readsnap(dir+file,'mass','dm')
    
    return mass


#def part2dens3d(part_pos, box_l, bin_x=128):
#    """
#    Calculate 3D matter density using numpy histograms
#    :param part_pos: particle positions in the shape (N, D), where N is particle number and D is dimension
#    :param box_l: box length in comoving Mpc/h
#    :param bin_x: desired bins per axis for the histogram
#    :return: density field
#    """
#    hist, _edges = np.histogramdd(np.vstack((part_pos[:, 0], part_pos[:, 1], part_pos[:, 2])).T,
#                                  bins=bin_x, range=[[0, box_l], [0, box_l], [0, box_l]])
#    del _edges
#    return hist
#
#def dens2overdens(density, mean_density=None):
#    """
#    Calculate the overdensity corresponding to a density field
#    :param density: input density field
#    :param mean_density: if defined normalisation is calculated according to (density - mean(density)) / mean_density
#    :return: overdensity field
#    """
#    assert np.ndim(density) == 3, 'density is not 3D'
#    if mean_density:
#        delta = (density - np.mean(density)) / mean_density
#    else:
#        delta = density / np.mean(density) - 1.
#    return delta
#
##
##def plot_ps_gadget(path_in,file_in,path_out,bin_x_,bin_k):
##    
##    print(path_in+'gadget_out/')
##    pos=ReadPos(path_in+'gadget_out/',file_in)
##    z=readheader(path_in+'gadget_out/'+file_in,'redshift')
##    
##
##    return
#    
#
#def power_spectrum(field_x, box_l, bin_k):
#    """
#        Measures the mass power spectrum of a 3D input field for a given number of bins in Fourier space.
#        :param field_x: 3D input field to compute the power spectrum of (typically the overdensity field), dimensionless
#        :param box_l: box length of image/cube/box or whatever, units of Mpc or Mpc/h
#        :param bin_k: number of bins in Fourier space
#        :return: power_k, k: 1D mass power spectrum of field_x, same units as [box_l]**3 and corresponding k values
#        """
#    assert np.ndim(field_x) == 3, 'field_x is not 3D'
#    box_pix = np.size(field_x, axis=0)  # pixel number per axis
#    box_dim = np.ndim(field_x)  # dimension
#
#    # This first 'paragraph' is to create masks of indices corresponding to one Fourier bin each
#    _freq = np.fft.fftfreq(n=box_pix, d=box_l/box_pix) * 2*np.pi
#    _rfreq = np.fft.rfftfreq(n=box_pix, d=box_l/box_pix) * 2*np.pi
#    _kx, _ky, _kz = np.meshgrid(_freq, _freq, _rfreq, indexing='ij')
#    del _freq, _rfreq
#    _k_abs = np.sqrt(_kx**2. + _ky**2. + _kz**2.)
#    del _kx, _ky, _kz
#    # The following complicated line is actually only creating a 1D array spanning k-space logarithmically from minimum _k_abs to maximum.
#    # To start slightly below the minimum and finish slightly above the maximum I use ceil and floor.
#    # To ceil and floor not to the next integer but to the next 15th digit, I multiply by 1e15 before flooring and divide afterwards.
#    # Since the ceiled/floored value is actually the exponent used for the logspace, going to the next integer would be way too much.
#    _k_log = np.logspace(np.floor(np.log10(np.min(_k_abs[1:]))*1.e15)/1.e15, np.ceil(np.log10(np.max(_k_abs[1:]))*1.e15)/1.e15, bin_k)
#
#    _field_k = np.fft.rfftn(np.fft.fftshift(field_x)) * (box_l/box_pix)**box_dim
#    power_k = np.empty(np.size(_k_log)-1)
#    for i in xrange(np.size(_k_log)-1):
#        mask = (_k_abs >= _k_log[i]) & (_k_abs < _k_log[i+1])
#        power_k[i] = np.mean(np.abs(_field_k[mask])**2.) / box_l**box_dim
#    del _k_abs, _field_k
#
#    k = (_k_log[1:] + _k_log[:-1]) / 2
#    del _k_log
#
#    return power_k, k

#def plot_ps_gadget(path_in,file_in,path_out,bin_x_,bin_k):
#    
##    print(path_in+'gadget_out/')
#    pos=ReadPos(path_in+'gadget_out/',file_in)
#    print(file_in)
#    z=readheader(path_in+'gadget_out/'+file_in,'redshift')
#    print(z)
#    L=readheader(path_in+'gadget_out/'+file_in,'boxsize')    
#    
##    pos=ReadPos('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/','out_z31p000')
##    z=readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/out_z31p000','redshift')
#    
#    density=part2dens3d(pos, box_l=L, bin_x=bin_x_)
#    delta = dens2overdens(density, np.mean(density))
#    P_k = power_spectrum(delta, np.float64(L), bin_k)[0]
#    k = power_spectrum(delta,  np.float64(L), bin_k)[1]
#            
#    plt.plot(k,P_k)
#    plt.xlabel(r"$k$ $[h \mathrm{Mpc}^{-1}]$")
#    plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^{3}]$")
#    plt.title('Gadget P(k), z=%2.f' %z)
#    plt.xscale('log')
#    plt.yscale('log')
#    #plt.show()
#    if not os.path.exists(os.path.dirname(path_out+'gadget/')):
#        os.makedirs(path_out+'gadget/')
#    plt.savefig(path_out+'gadget/ps_'+file_in+'.png', bbox_inches = "tight",dpi=300)
#    plt.close()
#    
#    return k,P_k
    
def plot_den_gadget(path_in,file_in,path_out,bin_x_,bin_k):
    
#    print(path_in+'gadget_out/')
    pos=ReadPos(path_in+'gadget_out/',file_in)
    vel=ReadVel(path_in+'gadget_out/',file_in)
#    mass=ReadMass(path_in+'gadget_out/',file_in)
    print(file_in)
    z=readheader(path_in+'gadget_out/'+file_in,'redshift')
    print(z)
    L=readheader(path_in+'gadget_out/'+file_in,'boxsize')    
    

    
#    pos=ReadPos('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/','out_z31p000')
#    z=readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/out_z31p000','redshift')
    
#    density=part2dens3d(pos, box_l=L, bin_x=bin_x_)
#    delta = dens2overdens(density, np.mean(density))
#    P_k = power_spectrum(delta, np.float64(L), bin_k)[0]
#    k = power_spectrum(delta,  np.float64(L), bin_k)[1]
#            
#    plt.plot(k,P_k)
#    plt.xlabel(r"$k$ $[h \mathrm{Mpc}^{-1}]$")
#    plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^{3}]$")
#    plt.title('Gadget P(k), z=%2.f' %z)
#    plt.xscale('log')
#    plt.yscale('log')
#    #plt.show()
#    if not os.path.exists(os.path.dirname(path_out+'gadget/')):
#        os.makedirs(path_out+'gadget/')
#    plt.savefig(path_out+'gadget/ps_'+file_in+'.png', bbox_inches = "tight",dpi=300)
#    plt.close()
    
    return pos,vel


#
#def ReadPos(dir,file):
#    pos=readsnap(dir+file,'pos','dm')
#    
#    return pos
#


file_in_picola,just_files_picola=pps.list_files_picola(sys.argv[1])   
file_in_gadget,just_files_gadget=pps.list_files_gadget(sys.argv[1])
file_in_ic,just_files_ic=pps.list_files_ic(sys.argv[1])

#
#for x in (just_files_picola):
#    print(x)
#    k_pic,P_k_pic=plot_ps_picola(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#    
#for x in (just_files_gadget):
#    print(x)    
#    k_gad,P_k_gad=plot_ps_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))


#print(just_files_ic)
#for x in just_files_ic:
#    k_ic,P_ic=plot_ps_ic(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#    plot_ps_ic_CAMB(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_ic,P_ic)
#    P_k2=plot_ps_ic_EH(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_ic,P_ic)
#
##for x, y in zip(just_files_gadget[::-1], just_files_picola):
#for x, y in zip(just_files_gadget, just_files_picola):    
#    
#    print(x, y)
#
##    print(x)    
#    k_gad,P_k_gad=plot_ps_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#    
##    print(y)
#    k_pic,P_k_pic=plot_ps_picola(sys.argv[1],y,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#  
##    print(readheader(sys.argv[1]+'gadget_out/'+x,'redshift'),readheader(sys.argv[1]+'picola_out/'+y,'redshift'))
#    plot_ps_comp(sys.argv[1],x,y,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_gad,P_k_gad,k_pic,P_k_pic)
##    
##   
##    print(x)
#    plot_ps_gadget_CAMB(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_gad,P_k_gad)
#    
##    print(x)
#    plot_ps_gadget_EH(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_gad,P_k_gad)
#    
#   
##    print(x,y)
#    plot_ps_picola_CAMB(sys.argv[1],y,x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
#    
##    print(x,y)
#    plot_ps_picola_EH(sys.argv[1],y,x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
#
##    print(x,y)
#    plot_ps_comp_all(sys.argv[1],x,y,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_gad,P_k_gad,k_pic,P_k_pic)
#
#
#print(just_files_gadget)
#
#for x in (just_files_gadget):
#    print(x)
#    k_gad,P_k_gad=plot_ps_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#    
##    k_pic,P_k_pic=plot_ps_picola(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
##    plot_ps_picola_CAMB(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
##    plot_ps_picola_EH(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
#
#




print(just_files_gadget)

for x in (just_files_gadget):
    print(x)
#    k_gad,P_k_gad=plot_ps_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
#    
    
    
    
    pos,vel=plot_den_gadget(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
    
    
    
    
    
    
##    k_pic,P_k_pic=plot_ps_picola(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]))
##    plot_ps_picola_CAMB(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
##    plot_ps_picola_EH(sys.argv[1],x,sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),k_pic,P_k_pic)
#


# import particle position and transfer to nbodykit using ArrayCatalog e to_mesh
    
#    to see the mesh insertion, see https://nbodykit.readthedocs.io/en/latest/mesh/creating.html#creating-mesh
    
    
    
    
    

pos
vel
#mass



#data_gad=pos;

data_test= []

# empty array
#arr = [] 

# some fake data
#data = numpy.ones(5, dtype=[
#      ('Position', ('f4', 3)),
#      ('Velocity', ('f4', 3))]
#      )

print(len(pos))




#read using 

print()


data = numpy.ones(262144, dtype=[
      ('Position', ('f4', 3)),
      ('Velocity', ('f4', 3))]
      )



data

#pos_to_nbk=data_test;
#
#pos_to_nbk
#
#data_test=array([ ['Position',pos],
#            ['Velocity',pos]]

# initialize a catalog directly from the structured array
#src = ArrayCatalog(data)

data_t= data.copy()
#
#data_pos=data['Position']
#
#data_pos

#data_t["Position"]=0

#data_t


data_t["Position"]=pos
data_t["Velocity"]=vel



data_t

##println()
#
#print("data_t shape:", data_t.shape)
#
#print("dtype:", data_t.dtype)




src_t= ArrayCatalog(data)

src_t
    
#src = ArrayCatalog(data)


src = ArrayCatalog(data_t)


## overwrite the Velocity column
#src['Velocity'] = src['Position'] + src['Velocity'] # 1 + 1 = 2
#
## overwrite the Position column
#src['Position'] = src['Position'] + src['Velocity'] # 1 + 2 = 3

print("Velocity = ", src.compute(src['Velocity'])) # all equal to 2
print("Position = ", src.compute(src['Position'])) # all equal to 3

#BoxSize=1.0
mesh = src.to_mesh(Nmesh=64,BoxSize=64)

#mesh = 


print("mesh = ", mesh)

plt.imshow(mesh.preview(axes=[0,1]))



#mesh stuff









# compute the power, specifying desired linear k-binning

#r = FFTPower(mesh, mode='1d', dk=0.005, kmin=0.01)
#r = FFTPower(mesh, mode='1d', dk=4, kmin=8)

r = FFTPower(mesh, mode='1d', dk=0.05, kmin=0.1)

#r

# the result is stored at "power" attribute
Pk = r.power
print(Pk)

Pk

print(Pk.coords)

Pkcoords=Pk.coords

Pkcoords



Pk_=Pk.coords

Pk_



#print(Pk.coords)

#print out the meta-data
for k in Pk.attrs:
    print("%s = %s" %(k, str(Pk.attrs[k])))
    
    
plt.figure()    
    
# print the shot noise subtracted P(k)
#plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'])

plt.loglog(Pk['k'], Pk['power'].real)
#plt.loglog(Pk['k'], Pk.attrs['shotnoise'])
Pk_k=Pk['k']

Pk_k

Pk_p=Pk['power']
Pk_p

Pk_pr=Pk['power'].real
Pk_pr

Pk_sn=Pk.attrs['shotnoise']

Pk_sn


#plt.imshow(Pk['k'], Pk['power'].real)

#plt.imshow(mesh.preview(axes=[0,1]))

#
#plt.imshow(mesh.preview(axes=[0,1]))
# format the axes
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
plt.xlim(0.1,2.5)









#from pmesh.pm import ParticleMesh, RealField, ComplexField


## a 8^3 mesh
#pm = ParticleMesh(Nmesh=[8,8,8])
#
## initialize a RealField
#rfield = RealField(pm)
#
## shape
#print("shape = ", rfield.shape)
#
#
#
#pm
#
#
## perform the forward FFT
#cfield = rfield.r2c()
#
#
## stores Nmesh/2+1 in z axis b/c of conjugate symmetry
#print("shape = ", cfield.shape)
#
## stores Nmesh/2+1 in z axis b/c of conjugate symmetry
#print("shape = ", cfield.shape)
#
## k=0 mode is the mean value of configuration space field
#print("mean of configuration space field from k=0 = ", cfield[0,0,0])
#
## perform the inverse FFT
#rfield2 = cfield.c2r()
#
## print the mean of the underlying array
#print("mean of real field = ", rfield2.value.mean())
#












#BoxSize=1.0

#print("mesh = ", mesh, "BoxSize = ",BoxSize)


#
#
#
## compute the power, specifying desired linear k-binning
#r = FFTPower(mesh, mode='1d', dk=0.005, kmin=0.01)
#
## the result is stored at "power" attribute
#Pk = r.power
#print(Pk)
#
#print(Pk.coords)
#
## print out the meta-data
#for k in Pk.attrs:
#    print("%s = %s" %(k, str(Pk.attrs[k])))
#    
#    
## print the shot noise subtracted P(k)
#plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'])
#
## format the axes
#plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
#plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
#plt.xlim(0.01, 0.6) 
#
## compute the 2D power
#r = FFTPower(mesh, mode='2d', dk=0.005, kmin=0.01, Nmu=5, los=[0,0,1])
#Pkmu = r.power
#print(Pkmu)
#
#print(Pkmu.coords)
#
#
## plot each mu bin
#for i in range(Pkmu.shape[1]):
#    Pk = Pkmu[:,i] # select the ith mu bin
#    label = r'$\mu$=%.1f' % (Pkmu.coords['mu'][i])
#    plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'], label=label)
#
## format the axes
#plt.legend(loc=0, ncol=2)
#plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
#plt.ylabel(r"$P(k, \mu)$ [$h^{-3}\mathrm{Mpc}^3$]")
#plt.xlim(0.01, 0.6)   
