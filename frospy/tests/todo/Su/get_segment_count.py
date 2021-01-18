# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 10:13:05 2018

@author: sujania
"""
import os
import sys
import glob
from frospy.core.segment import read

#if len(sys.argv) < 2:
#    print "USAGE python count_segment.py event.segments"
#    sys.exit()
#
#segment = sys.argv[1]
#seg_file=open(segment).read().splitlines()
#seg_count = len(seg_file)/4
#print 'Number of stations in segment file = %s'%seg_count
#print 'tw = %s'%seg_file[2].split()
#print 'fw = %s'%seg_file[1].split()

files = glob.glob("*.dat")
c = 0
for f in files:
    name = os.path.splitext(f)[0]
    seg = read(f)
    print ('Number of stations in %s is %s'%(name, len(seg)))
    c = c + len(seg)
print('total segments: %s'%c)
