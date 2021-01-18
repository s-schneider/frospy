#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:33:00 2019

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
from frospy.postprocessing.uncertainties import uncertainties_calculation
from frospy.postprocessing.AnalyseRuns.plot import misfits_cst, listdir
from frospy.postprocessing import plot
from frospy.postprocessing.plot import inv_summary, cst
from frospy.postprocessing.plot import plot_coeffs_per_mode
from frospy.postprocessing.read import read_inversion_summary
from frospy.postprocessing.misfit import plot_Qf
from frospy.core.splittingfunc import Set
#from frospy.core.splittingfunc import loadmodel # python2 branch
from frospy.splitting.load import loadmodel # python3 master branch
from frospy.core.splittingfunc.read import read_cst, get_cst, read_cst_S20RTS
from frospy.util.read import read_modes_in, get_mode_names
from frospy.util.base import sort_human 
from frospy.core.setup.settings import Setup, get_mode_list
from frospy.util.base import (chunking_list, split_digit_nondigit,
                            max_cc_degrees, max_sc_degrees)
from frospy.core.setup.settings import read as read_setup
from frospy.postprocessing.uncertainties import uncertainties_calculation
from frospy.core.database.write import write_cst
from frospy.core.database.query import get_db_tables, cst_query, db_query
from frospy.core.modes import read as read_modes


#modes_all = read_modes()
#mode = modes_all.select(name='0s0')[0]

##---------------------------------------------------------------------------

master_db_path = '/net/home/talavera/eejit/data/mycst.sqlite3'

host    = '/net/home/talavera/eejit'
cwd     = '/net/home/talavera/eejit/splitting'
#cwd     = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
#mode    = '02t07' #'03s00-08s02-09s02'#
#event ='cst+dst'#'%sevents'%'all'
#cpl   ='prem_TR' #R=0.2_AD
#mzero ='cst14+dst12' 
#name    = 'dst'

mode  = '10s02-11s02'
event ='%sevents'%'all'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='c00'
name    = 'QCC'
inv     = 18
damp    = 0.01

#event   = '%sevents'%'all'
#cpl     = 'cst+d00+cc' #cst+d00+cc #c00+d00 cst+d00
#mzero   = 'cst_02s00_c20=15 '#'08s02-09s02' #romanowicz begtramp
#mode  = '02s00-07s02'#
##mode  = '03s00-08s02-09s02'#
##event ='%sevents'%'good'
##cpl   ='cst+d00+cc' #c00+d00 cst+d00
##mzero ='08s02-09s02'
##name    = 'GC'
#inv     = 1
#damp    = 0.001

verbose = True
scut=4

cdir                = join(cwd,mode,event,cpl,mzero)
pickle_dir          = join(cdir, 'setup.json')
setup               = read_setup(pickle_dir)
setup.rundir        = join(cwd,mode,event,cpl,mzero)
db_path             = join(setup.rundir,'uncertainties','%s.sqlite3'%mode)

cerr = uncertainties_calculation(setup=setup, name=name, inv=inv,#host=host,
                               db=True, db_path=db_path, #model='data',
                               output_errors=False, jack_knive=True,
                               uncertainties_dir='uncertainties',
                               verbose=verbose,damp=damp,
                               )

cerr.plot(error_bars='plus',smax=scut)
s = Set()
s += cerr

#for m,coeffs in cerr.cst.items():
#    print(m)
#    for s,coe in coeffs.items(): 
#       print(m, s)
#       e_u =cerr.cst_errors[m][s]['upper_uncertainty']
#       e_l =cerr.cst_errors[m][s]['lower_uncertainty']
#       for c,l,u in zip(coe,e_l,e_u):
#           print(c,l,u)
           
#model     = None
#model     = ['S20RTS', 'AD']
#if model is not None:
#    for x in model:
#        s += loadmodel(format=x, setup=setup, 
#                                name = x,   
#        #                        name=x, 
#                                damp=0)  
#s.plot(error_bars='plus',smax=scut)

#write_cst(s, master_db_path, name, verbose=True, lcut=scut) 
#print get_db_tables(master_db_path)

SF = plot_coeffs_per_mode(SF_in=s, 
#                          label1 = "SC", 
                          mode='sc', 
                          plot_f_center=True, 
#                          plot_Q=True, kind = "dst", 
                          degree = 0, 
#                          ordering=["l", "sens"],
#                          l=0,
                          rot=0,
                          verbose=verbose,
                          )

SF = plot_coeffs_per_mode(SF_in=s, 
#                          label1 = "SC", 
                          mode='sc', 
#                          plot_f_center=True, 
                          plot_Q=True, kind = "dst", 
                          degree = 0, 
#                          ordering=["l", "sens"],
#                          l=0,
                          rot=0,
                          verbose=verbose,
                          )


from frospy.core.database.query import get_db_tables, cst_query, db_query
##db_path = '/net/home/talavera/codes/nmPy/nmpy/data/cst.sqlite3'
db_path = '/net/home/talavera/eejit/data/mycst.sqlite3'
model='QCC'
print(get_db_tables(db_path)) # Tables
c = cst_query(db_path, model, ['10S2'])
print(c)

#cerr.cst_errors["10S2"]["0"]['upper_uncertainty']=[10]
#cerr.cst_errors["10S2"]["0"]['lower_uncertainty']=[-10]
#cerr.cst_errors["10S2"]["0"]['uncertainty']=[10]
#cerr.cst_errors["11S2"]["0"]['upper_uncertainty']=[10]
#cerr.cst_errors["11S2"]["0"]['lower_uncertainty']=[-10]
#cerr.cst_errors["11S2"]["0"]['uncertainty']=[10]
#cerr.cst_errors["10S2-11S2"]["0"]['uncertainty']=[10]
#cerr.cst_errors["10S2-11S2"]["0"]['lower_uncertainty']=[-10]
#cerr.cst_errors["10S2-11S2"]["0"]['upper_uncertainty']=[10]