# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 17:15:56 2018

@author: sujania
"""

import os
import glob
import sys
from frospy.core.segment import read

seg = sys.argv[1]

files = glob.glob("*.%s"%seg)

for f in files:
    name = os.path.splitext(f)[0]
    seg = read(f)
    seg.sort()
    seg.write('%s.%s'%(name,seg))



#################################
    
import os
import glob
from frospy.core.segment import read

files = glob.glob("*.pickle")

for f in files:
    name = os.path.splitext(f)[0]
    seg = read(f)
    seg.sort()
    seg.write('%s.segmentZ'%(name))
