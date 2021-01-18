#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 10:09:16 2018

@author: talavera
"""
import glob
import numpy as np
from frospy.util.read import read_st
from frospy.core.segment import read
from os.path import basename, splitext

mode   = '03s00'
segs   = '%s-all' % mode
segdir = '//nfs/stig/talavera/radial-inversion/%s/segments/%s' % (mode, segs)

segments = glob.glob('%s/[0-9]*[A-Z].dat' % (segdir))
events   = [splitext(basename(seg))[0] for seg in segments]

for sfile, ev in zip(segments, events):
    ah = '//nfs/stig/deuss/modes/alldatafiles/%s.ahx' % ev
    st = read_st(ah, 'ah')
    
    seg = read(sfile)
    
    for s in seg:
        try:
            st.select(station = s.station)[0]
        except Exception:
            print 'error:', ev, s.station