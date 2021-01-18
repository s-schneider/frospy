#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 13:39:56 2018

@author: talavera
"""

from frospy.postprocessing.selection import outlier_removal

#mdir   = '/Users/sujania/Documents/PhD/UU.stig/radial-inversion'
#mdir  = '/net/home/talavera/eejit/splitting'
mdir  = '/scratch/talavera/splitting'
#mode  = '03s00-08s02-09s02'
#event ='%sevents'%'all'
#cpl   ='cst+d00+cc' #c00+d00 cst+d00
#mzero ='08s02-09s02' #romanowicz cst+d00_best_result
mode  = '00s15'
event ='cst'#'%sevents'%'all'
cpl   ='prem_TRZ' #c00+d00 cst+d00
mzero ='cst20' #romanowicz cst+d00_best_result
run_dir = os.path.join(mdir, mode, event, cpl, mzero)
db_path = os.path.join(run_dir,'inversion_out/inversion.db')

damp = 0.01
iteration = 10
max_misfit = 1
array=True
print_seg=False
seg_suffix = 'segment'
verbose=True
#sets=["VHZ"]#'all'
sets='all'
outlier_removal(run_dir=run_dir, damp=damp, iteration=iteration, 
                max_misfit=max_misfit, array=array, print_seg=print_seg, 
                db_path=db_path, seg_suffix=seg_suffix, sets=sets,
                verbose=verbose)

#import os
#import sys
#import glob
#import numpy as np
#from os.path import join, basename, splitext
#from frospy.core.segment import read
#from obspy.core import AttribDict
#from frospy.core.database.query import db_query
#from frospy.core.database.query import get_db_tables
#from frospy.core.segment import read
#
#
#def outlier_removal(rundir, db_path, damp, iteration=10, max_misfit=1.0,
#                    array=False, print_seg=False, seg_suffix = 'dat',
#                    overwrite=False):
#    """
#    array: PBS array style of inversion
#    print_seg: print new segments files
#
#    maindir = '//nfs/stig/talavera/radial-inversion/02s00/allevents/'
#    rundir = '%s/cst+d00/cst_02s00/rerun/' % maindir
#    db_path ='%s/inversion_out/inversion.db' % rundir
#    damp = 0.001
#    iteration = 10
#    max_misfit = 1.0
#    array=True
#    print_seg=False
#    seg_suffix = 'dat'
#    outlier_removal(rundir, db_path, damp, iteration, max_misfit,
#                    array, print_seg, seg_suffix)
#    """
#
#    # Reading the segments names
#    tables = get_db_tables(db_path)
#    sets_name = []
#    for t in tables:
#        try:
#            n = float(t)
#        except:
#            sets_name.append(t)
#    sets_name.remove('initial')
#
#    # total number of segments counter
#    old = 0
#    new = 0
#    new_set = 0
#    set_misfit = 0
#    
#    if len(sets_name) == 1:
#        allevents = set(db_query(db_path, sets_name[0],
#                        condition_dict={'iter': 1, 'damp': damp},
#                        select='event'))
#        allevents = [e[0] for e in allevents]
#        for event in allevents:
#            # reading segments
#            if array == False:
#                seg = read('%s/%s.%s' % (rundir, event, seg_suffix))
#            else:
#                seg = read('%s/%s/%s.%s' % (rundir, event, event, seg_suffix))
#
#            # reading stations misfits from inversion output
#            stations = db_query(db_path, sets_name[0],
#                                condition_dict={'event': event,
#                                                'iter': iteration,
#                                                'damp': damp},
#                                select='station, misfit')
#
#            old  += len(stations)
#            rm_sta = []
#            # removing stations with a high misfit
#            for sta, mf in stations:
#                set_misfit += mf
#                if mf > max_misfit:
#                    rm_sta.append(sta)
#            for sta in rm_sta:
#                seg.remove(seg.select(station=sta)[0])
#            new += len(seg)
#
#            # printing segments
#            if print_seg == True:
#                if array == False:
#                    seg.write('%s/%s.%s' % (rundir, event, seg_suffix), 
#                              overwrite=overwrite)
#                else:
#                    seg.write('%s/%s/%s.%s' % (rundir, event,
#                                               event, seg_suffix), 
#                              overwrite=overwrite)
#
#        set_misfit = set_misfit/ new
#        print(sets_name[0], set_misfit)         
# 
#    else:
#        for sset in sets_name:
#            allevents = set(db_query(db_path, sset,
#                            condition_dict={'iter': 1, 'damp': damp},
#                            select='event'))
#            allevents = [e[0] for e in allevents]
#            for event in allevents:
#                # reading segments
#                seg = read('%s/%s/%s/%s.%s' % (rundir, sset, event,
#                                               event, seg_suffix))
#
#                # reading stations misfits from inversion output
#                stations = db_query(db_path, sset,
#                                    condition_dict={'event': event,
#                                                    'iter': iteration,
#                                                    'damp': damp},
#                                    select='station, misfit')
#
#                old  += len(stations)
#                rm_sta = []
#                # removing stations with a high misfit
#                for sta, mf in stations:
#                    set_misfit = set_misfit + mf
#                    if mf > max_misfit:
#                        rm_sta.append(sta)
#                for sta in rm_sta:
#                    seg.remove(seg.select(station=sta)[0])
#                new += len(seg)
#                new_set += len(seg)
#
#                # printing segments
#                if print_seg == True:
#                    seg.write('%s/%s/%s/%s.%s' % (rundir, sset, event,
#                                                  event, seg_suffix), 
#                              overwrite=overwrite)
#            print(sset, set_misfit,new_set)
#            set_misfit = set_misfit/ new_set
#            print(sset, set_misfit) 
#            set_misfit = 0
#            new_set = 0
#
#    print 'Original segments: %s' % (old)
#    print '     New segments: %s, with misfit below %g' % (new, max_misfit)
#    return

