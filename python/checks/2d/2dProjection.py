#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:44:01 2023

@author: asus
"""

def plot_2d_proj(mesh):
    
    from matplotlib import pyplot as plt
    
    # plt.figure()    
    plt.imshow(mesh.preview(axes=[0,1]))


    