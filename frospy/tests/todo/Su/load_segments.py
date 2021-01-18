#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:18:53 2019

@author: talavera
"""
import os
from os.path import join
from frospy.core.segment import Segment
from frospy.core.segment import read as read_segment

# segments need to be upper case
#n_all=['00','00','00','00','01','01','01','01','01','02',
#       '02','02','02','03','03','04','05','05','05','05',
#       '06','07','07','08','08','09','11','12','12','12']
#
#l_all=['02','05','06','07','04','07','08','09','10','01',
#       '06','12','13','06','09','09','03','07','08','12',
#       '10','05','07','06','07','11','09','06','07','13']
#
#n_all=['02','02','05','04',
#       '06','07','08','08','09','11','12','12','12']
#
#l_all=['12','13','12','09',
#       '10','07','06','07','11','09','06','07','13']

min_snr=2
n_all=['00']
l_all=['24']
for (n,l) in zip(n_all,l_all):
    mode = '%sS%s'%(int(n),int(l))
    cwd = '/net/home/talavera/eejit/splitting/%ss%s/segments/'%(n,l)
    print(mode)
    # Z component
    comp  = 'Z'
    seg_Z = read_segment(ifile='db', modes=mode, channel="VH%s"%comp, min_snr=min_snr)
    allevents = []
    for seg in seg_Z:
        if seg.stats.event not in allevents:
            allevents.append(seg.stats.event)
    print(mode, comp, len(allevents), len(seg_Z))
    
    for e in allevents:
        if e[-1] != 'Z' and e[-1] != 'Q':
            s = Segment()
            for seg in seg_Z:
                if seg.stats.event == e:
                    s.append(seg)
            print e, len(s)
        
            cdir = join(cwd, 'segments_%s'%comp)
            s.sort()
            s.write('%s/%s.segment%s'%(cdir,e,comp))
    
    
    # T component
    comp  = 'T'
    seg_T = read_segment(ifile='db', modes=mode, channel="VH%s"%comp, min_snr=min_snr)
    allevents = []
    for seg in seg_T:
        if seg.stats.event not in allevents:
            allevents.append(seg.stats.event)
    print(mode, comp, len(allevents), len(seg_T))
    
    for e in allevents:
        if e[-1] != 'Z' and e[-1] != 'Q':
            s = Segment()
            for seg in seg_T:
                if seg.stats.event == e:
                    s.append(seg)
            print e, len(s)
        
            cdir = join(cwd, 'segments_%s'%comp)
            s.sort()
            s.write('%s/%s.segment%s'%(cdir,e,comp))
    
    
    # R component
    comp  = 'R'
    seg_R = read_segment(ifile='db', modes=mode, channel="VH%s"%comp, min_snr=min_snr)
    allevents = []
    for seg in seg_R:
        if seg.stats.event not in allevents:
            allevents.append(seg.stats.event)
    print(mode, comp, len(allevents), len(seg_R))
    
    for e in allevents:
        if e[-1] != 'Z' and e[-1] != 'Q':
            s = Segment()
            for seg in seg_R:
                if seg.stats.event == e:
                    s.append(seg)
            print e, len(s)
        
            cdir = join(cwd, 'segments_%s'%comp)
            s.sort()
            s.write('%s/%s.segment%s'%(cdir,e,comp))
