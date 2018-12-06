# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:47:51 2017
@author: eloisechakour
"""
import numpy
from pygadgetreader import *

readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/snapshot_000','time',debug=1)
a=readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/snapshot_000','redshift')
print(a)
#def test():
#
#    readheader('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/snapshot_000','redshift')
#
#    print(a)
#    return a