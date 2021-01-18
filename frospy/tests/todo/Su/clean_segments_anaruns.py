#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 13:30:42 2018

@author: talavera
"""
import os
import sys
import numpy as np

#if len(sys.argv) < 3:
#    print "USAGE python clean_segments.py it_max max_mf"
#    sys.exit()
#
#it_max = sys.argv[1]
#max_mf = sys.argv[2]

it_max = 9
max_mf = 1
cwd            = os.getcwd()
config_folder  = '%s/config'%cwd
seg_folder_old = '%s/config/segments_Z.highmf'%cwd
seg_folder_new = '%s/config/segments_Z'%cwd
data_folder    = '/home/talavera/broadband-syn/0-7.5mHz'
#data_folder     = '%s/config/ah_dir'%(cwd)
syn_folder     = '%s/Iter_%s/Z/synseis'%(cwd, it_max)
events_file    = '%s/config/events_list.dat'%cwd

events  = open(events_file).read().splitlines()[1::]
seg_old = 0
seg_new = 0

for e in events:  
    data_file = '%s/%s.ahx'%(data_folder, e)
    syn_file  = '%s/%s.ahx.syn'%(syn_folder, e)
    seg_file  = '%s/%s.dat'%(seg_folder_old, e)
    
    seg = spectrum(data=data_file, syn=syn_file, segment_file=seg_file, 
                   runall=True, max_mf=max_mf)
    seg.write('%s/%s.dat'%(seg_folder_new, e))
    
    seg_old = len(open(seg_file).read().splitlines())/4 + seg_old
    seg_new = len(seg) + seg_new

VR = open('%s/Total_VR.out'%syn_folder).read().splitlines()[1].split(':')[1]
mf = 1 - float(VR)/100    
print 'Original  segment count: %s, Old misfit: %s'%(seg_old,mf)
print '     New  segment count: %s, with misfit below %s'%(seg_new, max_mf)
    
