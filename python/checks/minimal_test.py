#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 19:03:45 2023

@author: asus
"""



# import nbodykit

import numpy as np

import matplotlib
matplotlib.use('agg')   #deal with figures wihotut winow forwards (non iteractive, like slurm)
from matplotlib import pyplot as plt


# plt.figure()  
k = np.random.random ([3,4]) * 10
plt.imshow(k)
plt.savefig("/home/cunhad/project/cunhad/production/NBodyWake/python/checks/test.png", bbox_inches = "tight",dpi=300)
plt.close()

