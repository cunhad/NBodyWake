
# coding: utf-8

# In[1]:

#matplotlib inline


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

#uniform particle distribution

#data_path_n="/home/acer/Documents/storage/irulan/test/uniform/50Mpc_128c_zi65_nowakes/Analysis/stat/box_statistics/healpy/*64_c.txt"


# In[7]:

#uniform particle distribution with pole

#data_path_n="/home/acer/Documents/storage/irulan/test/uniform/50Mpc_128c_zi65_nowakes/Analysis/stat/box_statistics/healpy/*8_pole.txt"


# In[28]:

#from laptop
#data_path="/home/acer/Documents/storage/40Mpc_192c_zi65_wakeGmu1t10m5zi40s/Analysis/stat/box_statistics/healpy/*8_pole.txt"
#data_path_n="/home/acer/Documents/storage/40Mpc_192c_zi65_nowakes/Analysis/stat/box_statistics/healpy/*8_pole.txt"


# In[59]:

#from laptop wake centered
#data_path="/home/acer/Documents/storage/40Mpc_192c_zi65_wakeGmu1t10m5zi40s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path_n="/home/acer/Documents/storage/40Mpc_192c_zi65_nowakes/Analysis/stat/box_statistics/healpy/*64_c.txt"


# In[24]:

#from irulan 32Mpc_512c_zi63
#data_path="/home/acer/Documents/storage/irulan/32Mpc_512c_zi63_wakeGmu1t10m6zi15s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path2="/home/acer/Documents/storage/irulan/32Mpc_512c_zi63_wakeGmu1t10m6zi7s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path_n="/home/acer/Documents/storage/irulan/32Mpc_512c_zi63_nowakes/Analysis/stat/box_statistics/healpy/*64_c.txt"


# In[9]:

#from irulan 50Mpc_512c_zi65 interesting

#data_path="/home/acer/Documents/storage/irulan/interesting/50Mpc_512c_zi65_wakeGmu6t10m6zi40s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path2="/home/acer/Documents/storage/irulan/interesting/50Mpc_512c_zi65_wakeGmu6t10m6zi20s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path3="/home/acer/Documents/storage/irulan/interesting/50Mpc_512c_zi65_wakeGmu3t10m6zi40s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path_n="/home/acer/Documents/storage/irulan/interesting/50Mpc_512c_zi65_nowakes/Analysis/stat/box_statistics/healpy/*64_c.txt"


# In[16]:

#from irulan 16Mpc_512c_zi63

#data_path="/home/acer/Documents/storage/irulan/16Mpc_512c_zi63_wakeGmu1t10m6zi15s/Analysis/stat/box_statistics/healpy/*64_c.txt"
#data_path_n="/home/acer/Documents/storage/irulan/16Mpc_512c_zi63_nowakes/Analysis/stat/box_statistics/healpy/*64_c.txt"



# In[9]:

#from guillimim 8 nodes simplest one 32Mpc_96c_zi63_nowakes

data_path="/home/asus/Dropbox/extras/storage/guillimin/old/32Mpc_96c_zi63_nowakes/Analysis/stat/box_statistics/max_1dproj/*64.txt"


# In[13]:

#import glob
#files_list_n = sorted(glob.glob(data_path_n))
#files_list_n


# In[33]:

#max_density=np.loadtxt(files_list_n[0], delimiter=' ');
#max_density
#hp.mollview(max_density,title='50Mpc_512c_zi65_nowakes_z2')


# In[10]:

import glob
files_list = sorted(glob.glob(data_path))
files_list


# In[21]:

max_density=np.loadtxt(files_list[1], delimiter=' ');
#max_density.0p
#hp.mollview(max_density,title='50Mpc_512c_zi65_wakeGmu6t10m6zi40s_z2')
plot = hp.mollview(max_density)


# In[30]:

from pylab import *
import matplotlib.pyplot as plt


# In[44]:

hp.mollview(max_density)
plt.savefig('test3.pdf')


# In[17]:

#import glob
#files_list2 = sorted(glob.glob(data_path2))
#files_list2


# In[35]:

#max_density=np.loadtxt(files_list2[0], delimiter=' ');
#max_density.0p
#hp.mollview(max_density,title='50Mpc_512c_zi65_wakeGmu6t10m6zi20s_z2')


# In[19]:

import glob
#files_list3 = sorted(glob.glob(data_path3))
#files_list3


# In[36]:

#max_density=np.loadtxt(files_list3[0], delimiter=' ');
#max_density.0p
#hp.mollview(max_density,title='50Mpc_512c_zi65_wakeGmu3t10m6zi40s_z2')


# In[ ]:



