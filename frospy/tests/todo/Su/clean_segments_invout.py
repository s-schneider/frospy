s#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:05:11 2018

@author: talavera
"""
import os
import sys
import glob
import numpy as np
from os.path import join, basename, splitext
from frospy.core.segment import read
from obspy.core import AttribDict

    
def read_dot_misfit(file):
    """
    Reads all cmt.misfit files, such as 060994A.misfit, found in folder
    and returns a list of the misfits

    returns:
    list: seg_misfits
    """
    # Calculating segment misfits
    event = AttribDict()
    cmt = file.split('/')[-1].split('.')[0].split('_')[0]
    staname = []
    misfit = []
    with open(file, 'r') as fh:
        content = fh.readlines()
    for line in content:
        if len(line.split()) == 1:
            staname += [str(line.split()[0])]
        if len(line.split()) == 3:
            misfit += [float(line.split()[-1])]
    stamf = AttribDict()
    for sta, mf in zip(staname, misfit):
        stamf[sta] = mf
    event[cmt] = stamf
    return event

#mdir   = '/Users/sujania/Documents/PhD/UU.stig/radial-inversion'
mdir  = '//nfs/stig/talavera/radial-inversion'
mode  = '05s00'
event = '%sevents'%'deep'
cpl   = 'cst+d00+cc' 
mzero = 'cst_05s00_13s02' 
damp  = '%4.2e' % (0.0001) # input

run_dir  = os.path.join(mdir, mode, event, cpl, mzero, 'rerun') # input

#seg_sets = ['inversion_out'] # input
seg_sets = ['05s00','13s02'] # input
seg_mif  = []
max_mf   = 0.5 # input
seg_dir = []
segments = []
for sset in seg_sets:
    seg_mif += glob.glob('%s/%s/[0-9]*[A-Z]*d%s.misfit' % (run_dir, 
                                                            sset, damp))
    seg_dir.append(os.path.join(run_dir,sset)) # input
    
#seg_dir  = run_dir # input
for sset in seg_dir:
    segments.extend(glob.glob('%s/[0-9]*[A-Z].dat' % (sset)))

segments.sort()
seg_mif.sort()
old = 0
new = 0
for sfile, ffile in zip(segments, seg_mif):
    seg  = read(sfile)
    old  += len(seg)
    misfit = read_dot_misfit(ffile)
    rm_sta = []
    for event, stations in misfit.iteritems():
        for sta, mf in stations.iteritems():
            if mf > max_mf:
                rm_sta.append(sta)
    for sta in rm_sta:
        seg.remove(seg.select(station=sta)[0])
    new += len(seg)
#    seg.write('%s/%s.dat' % (seg_dir, event))

print 'Original  segment count: %s' % (old)
print '     New  segment count: %s, with misfit below %s' % (new, max_mf)
    
    

