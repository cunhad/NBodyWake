
# coding: utf-8

# In[1]:

import numpy as np
import healpy as hp


# In[2]:

NSIDE = 4096
NPix=hp.nside2npix(NSIDE)
NPix


# In[3]:

angles=hp.pix2ang(NSIDE,range(NPix))
#angles


# In[10]:

np.savetxt('~/Documents/angles4096_t.txt',angles[0], delimiter=' ')
np.savetxt('~/Documents/angles4096_p.txt',angles[1], delimiter=' ')


# In[5]:

angles[1]


# In[6]:

#angles.tofile('angles128_test.cvs')


# In[ ]:


#sampl = np.random.uniform(low=-1.0, high=1.0, size=(NPix,))


# In[ ]:

#np.savetxt('random_test8.txt',sampl, delimiter=' ')


# In[ ]:



