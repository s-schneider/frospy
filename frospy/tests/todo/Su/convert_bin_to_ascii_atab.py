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

HEADER_FILE = 'atab.header'
ATA_FILE    = 'atab.mat'
RHS_FILE    = 'atab.rhs'

###################
# READ ATAB FILES #
###################

ata = np.fromfile(ATA_FILE, dtype=np.dtype('f8')) 
rhs = np.fromfile(RHS_FILE, dtype=np.dtype('f8'))

with open(HEADER_FILE, mode='rb') as fhead: # b is important -> binary
    fileContent = fhead.read()

M = struct.unpack('i', fileContent[0:4])[0] # model parameter size
misfit = struct.unpack('d', fileContent[4:12])[0]

ata = ata.reshape([M,M])
rhs = rhs.reshape([M,1])

np.savetxt('ascii_atab.mat',ata,fmt='%16.4e')#,delimiter='\n')
np.savetxt('ascii_atab.rhs',rhs,fmt='%16.4e',delimiter='\n')

print "%d %.4e" %(M,misfit)
