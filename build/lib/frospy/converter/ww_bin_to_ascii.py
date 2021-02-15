#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:04:29 2017

@author: karaoglu
"""

import numpy as np
import sys

if len(sys.argv) < 2:
    print("USAGE python bin_to_ascii.py binary_file")
    sys.exit()

INFILE = sys.argv[1]
###################
# READ BIN FILES #
###################

mat = np.fromfile(INFILE, dtype=np.dtype('f8'))
np.savetxt('%s.ascii' % INFILE, mat, fmt='%16.4e', delimiter='\n')
