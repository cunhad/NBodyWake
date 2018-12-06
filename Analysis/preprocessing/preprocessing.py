#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 12:09:26 2018

@author: Disrael Cunha
"""


#path_in='/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/'
#file_name='snapshot'

import numpy as np
#import sys
#import natsort
import glob

import re 

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#file_in=glob.glob(sys.argv[1]+sys.argv[2]+'*')

#for x in sorted_nicely(file_in):
#    print(x)


def list_files_gadget(path_files):
    file_in=glob.glob(path_files+'gadget_out/snapshot*')
    file_in=sorted_nicely( file_in )
    just_files=[x.split('/')[-1] for x in file_in]
    return file_in,just_files

def list_files_picola(path_files):
    file_in=glob.glob(path_files+'picola_out/out*.0')
    file_in=[x[:-2] for x in file_in]
    file_in=sorted_nicely( file_in )
    file_in=file_in[::-1]
    just_files=[x.split('/')[-1] for x in file_in]
    return file_in,just_files

##file_in,just_files=list_files_gadget('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/')
#file_in,just_files=list_files_picola('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/')
#for x in (just_files):
#    print(x)