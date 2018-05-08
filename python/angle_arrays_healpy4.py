
# coding: utf-8

# In[1]:

matplotlib inline


# In[2]:

import numpy as np
import healpy as hp


# In[3]:

angles=np.loadtxt('angles64.txt', delimiter=' ');
#angles


# In[4]:

NPix=angles.shape[1]
NSIDE = np.int(np.sqrt(NPix/12))


# In[5]:

NSIDE,NPix


# In[6]:

#from guillimim 8 nodes simplest one 32Mpc_96c_zi63_nowakes

data_path="/home/asus/Dropbox/extras/storage/guillimin/old/32Mpc_96c_zi63_nowakes/Analysis/stat/box_statistics/max_1dproj/*64.txt"


# In[7]:

import glob
files_list = sorted(glob.glob(data_path))
files_list


# In[9]:

max_density=np.loadtxt(files_list[1], delimiter=' ');
#max_density.0p
#hp.mollview(max_density,title='50Mpc_502c_zi65_wakeGmu6t10m6zi40s_z2')
plot = hp.mollview(max_density)


# In[10]:

from pylab import *
import matplotlib.pyplot as plt


# In[11]:

hp.mollview(max_density)
plt.savefig('test2.pdf')


# In[ ]:



