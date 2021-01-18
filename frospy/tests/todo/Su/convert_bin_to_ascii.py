#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:04:29 2017

@author: karaoglu
"""

import numpy as np
import struct
import sys
import os
import re

if len(sys.argv) < 2:
    print "USAGE python bin_to_ascii.py binary_file M(model parameter space length)"
    sys.exit()

INFILE = sys.argv[1]
M      = int(sys.argv[2])

###################
# READ ATAB FILES #
###################

mat = np.fromfile(INFILE, dtype=np.dtype('f8')) 

N = np.sqrt(len(mat))

if N != M:
    print "Expected: %d, Received: %d" %(N,M)
    sys.exit()
mat = mat.reshape([M,M])

np.savetxt('%s.ascii'%INFILE,mat,fmt='%16.4e',delimiter='\n')
