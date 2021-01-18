#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:54:43 2019

@author: talavera
"""
import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.postprocessing.AnalyseRuns.plot import misfits_cst, listdir
from frospy.postprocessing.misfit import plot_Qf
from frospy.core.splittingfunc import Set
from frospy.core.splittingfunc import loadmodel # python2 branch
#from frospy.splitting.load import loadmodel # python3 master branch
from frospy.core.splittingfunc.read import read_cst, get_cst, read_cst_S20RTS
from frospy.util.read import read_modes_in, get_mode_names
from frospy.util.base import sort_human 
from frospy.core.setup.settings import Setup, get_mode_list
from frospy.util.base import (chunking_list, split_digit_nondigit,
                            max_cc_degrees, max_sc_degrees)
from frospy.core.database.write import write_cst
#
model = 'GD'
cdir = '/net/home/talavera/eejit/splitting/'
#ifiles = glob.glob(join(cdir,'**s00/measurements','%s.*'%model))
#ifiles.append(glob.glob(join(cdir,'02s00-07s02/measurements','%s.*'%model))[0])
ifiles= glob.glob(glob.glob(join(cdir,'11s00-27s02/measurements','%s.*'%model))[0])
sort_human(ifiles)
sf = Set()

for f in ifiles:  
    mode = f.split('/')[6]
    input = OrderedDict() #AttribDict()
    input['rundir'] = f
    input['startmodel'] = 'PREM'
    if '-' not in mode:
        print mode
        input['modes'] = OrderedDict([(mode, 0)])
        setup=Setup('CST', input)
        sf += loadmodel(setup=setup, ifile = f, 
                        name = model, damp = 0)
    else:
        mode1 = mode.split('-')[0]
        mode2 = mode.split('-')[1]
        print mode1, mode2
        input['modes'] = OrderedDict([(mode1, 0), (mode2, 0)])        
        setup=Setup('CST', input)
        sf += loadmodel(setup=setup, ifile = f, 
                        name = model, damp = 0)


#model = 'REM'
#cdir = '/net/home/talavera/eejit/data/%s/*.cst'%model
#ifiles = glob.glob(cdir)
#sort_human(ifiles)
#sf = Set()
#
#for f in ifiles:  
#    mode = f.split('/')[-1].split('.')[0]
#    input = OrderedDict() #AttribDict()
#    input['rundir'] = f
#    input['startmodel'] = 'PREM'
#    print mode
#    input['modes'] = OrderedDict([(mode, 2)])
#    setup=Setup('CST', input)
#    sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0)
#        
 
#model = '1.20'
#cdir = '/net/home/talavera/eejit/data/lythgoe/%s/*.cst'%model
#ifiles = glob.glob(cdir)
#sort_human(ifiles)
#sf = Set()
#model="$\phi$=%s"%model
#for f in ifiles:  
#    mode = f.split('/')[-1].split('.')[0]
#    input = OrderedDict() #AttribDict()
#    input['rundir'] = f
#    input['startmodel'] = 'PREM'
#    print mode
#    input['modes'] = OrderedDict([(mode, 0)])
#    setup=Setup('CST', input)
#    sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0)
        
 
   
# sf = Set()
# ic = 'begtramp'; IC='BT'
# #mtype= 'S20RTS+CRUST+'
# mtype= ''
# model = '%s%s'%(mtype,IC)
# radial = ['02s00','05s00','11s00']
# innerc = ['07s02','13s02','27s02']
# cdir1 = '/net/home/talavera/radial-inversion'
# cdir2 = 'synthetics/IC/mcst-%s%s.dat'%(mtype,ic)
# cst_files = [join(cdir1,m,cdir2) for m in radial]

# sf = Set()
# for i, f in enumerate(cst_files):
#     input = OrderedDict() #AttribDict()
#     input['rundir'] = cst_files[i]
#     input['modes'] = OrderedDict([(radial[i], 0), (innerc[i], 4)])
#     input['modes_cc'] = OrderedDict([('%s-%s'%(radial[i], innerc[i]),[2, 2])])
#     input['startmodel'] = 'PREM'
#     setup=Setup('CST', input)
#     sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0)    

# radial = ['03s00','04s00']
# innerc1 = ['08s02','10s02']
# innerc2 = ['09s02','11s02']
# cdir1 = '/net/home/talavera/radial-inversion'
# cdir2 = 'synthetics/IC/mcst-%s%s.dat'%(mtype,ic)
# cst_files = [join(cdir1,m,cdir2) for m in radial]

# for i, f in enumerate(cst_files):
#     input = OrderedDict()#AttribDict()
#     input['rundir'] = cst_files[i]
#     input['modes'] = OrderedDict([(radial[i], 0), 
#                                   (innerc1[i], 4), 
#                                   (innerc2[i], 4)])
#     input['modes_cc'] = OrderedDict([('%s-%s'%(radial[i], innerc1[i]),[2, 2]),
#                                      ('%s-%s'%(radial[i], innerc2[i]),[2, 2]),
#                                      ('%s-%s'%(innerc1[i],innerc2[i]),[0, 4])])
#     input['startmodel'] = 'PREM'
#     setup=Setup('CST', input)
#     sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0)  



# radial  = ['03s00','04s00']
# innerc1 = ['08s02','10s02']
# innerc2 = ['09s02','11s02']
# modes   = zip(radial,innerc1,innerc2)
# c1 = '/net/home/talavera/radial-inversion'
# c2 = 'synthetics/IC/'
# c3 = 'mcst-%s%s.dat'%(mtype,ic)

# cst_files = [join(c1,m1,c2,"%s-%s"%(m2,m3),c3) for m1,m2,m3 in modes]

# for i, f in enumerate(cst_files):
#     input = OrderedDict()#AttribDict()
#     input['rundir'] = cst_files[i]
#     input['modes'] = OrderedDict([(innerc1[i], 4), 
#                                   (innerc2[i], 4)])
#     input['modes_cc'] = OrderedDict([('%s-%s'%(innerc1[i],innerc2[i]),[0, 4])])
#     input['startmodel'] = 'PREM'
#     setup=Setup('CST', input)
#     sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0) 
  


#sf = Set()
#model = 'GLW'
#radial  = ['04s00']
#innerc1 = ['10s02']
#innerc2 = ['11s02']
#modes   = zip(radial,innerc1,innerc2)
#c1 = '/net/home/talavera/radial-inversion'
#c2 = 'synthetics/IC/'
#c3 = 'mcst-%s.dat'%model
#
#cst_files = [join(c1,m1,c2,"%s-%s"%(m2,m3),c3) for m1,m2,m3 in modes]
#
#for i, f in enumerate(cst_files):
#    input = OrderedDict()#AttribDict()
#    input['rundir'] = cst_files[i]
#    input['modes'] = OrderedDict([(innerc1[i], 4)])
#    input['startmodel'] = 'PREM'
#    setup=Setup('CST', input)
#    sf += loadmodel(setup=setup, ifile = f, name = model, damp = 0) 
#  
db_path = '/net/home/talavera/codes/nmPy/nmpy/data/cst.sqlite3'
write_cst(sf, db_path, model, verbose=True)    

from frospy.core.database.query import get_db_tables, cst_query, db_query
print get_db_tables(db_path)
#c = cst_query(db_path, model, ['13S2'])
#c = cst_query(db_path, model, ['7S2'])
c = cst_query(db_path, model, ['3S0'])
print c
#for ce in c:
#   print ce[0], ce[1], ce[2], float(ce[3].split()[0]), float(ce[3].split()[1])