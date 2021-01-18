#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:56:02 2019

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
from frospy.postprocessing.plot import plot_coeffs_per_mode

# 2 modes
cwd   = '/net/home/talavera/eejit/splitting'
#cwd   = '/Users/sujania/Documents/PhD/UU.eejit/splitting'

#mode = '01s00-04s02-00s10'#
#mode  = '02s00-07s02'
#mode  = '03s00-08s02-09s02'
#mode  = '08s02-09s02'
#mode  = '04s00-10s02-11s02'
mode  = '10s02-11s02'
#mode  = '05s00-13s02'
#mode  = '11s00-27s02'

event ='%sevents'%'all'
cpl   ='cstC+dst' #c00+d00 cst+d00 cst+d00+cc
#mzero ='08s02-09s02' # begtramp cst+d00
#mzero ='10s02-11s02' # begtramp cst+d00
#mzero ='good' # begtramp cst+d00
mzero="c00"
inv   = 1#
enditer= 25#setup.iterations[1]
damp   = 0.01

pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.json')
setup     = read_setup(pickle_dir)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
rundir    = join(cwd,mode,event,cpl,mzero)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
dampings  = setup.damping
startiter = setup.iterations[0]
iters     = range(startiter, enditer+1)
#model     = ['S20RTS', 'AD']
model     = ['S20RTS']

##### best result
#sf = Set()
#if inv == 1:
#    invs = 'inversion_out'
#else:
#    invs = 'inversion_out_%s'%inv
#    
#cst_files = glob.glob('%s/%s/mnew-it%s-d%.2e.dat' % 
#                      (rundir, invs, enditer, damp))
#for i, f in enumerate(cst_files):
#    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
#    sf += loadmodel(setup=setup, ifile = f, 
#                    name ="Our inversion", 
##                    name ="it=%s, d=%s"%(str(it),damp), 
#                    damp = damp)
#sf.plot_map(cmap="su", lon_0=-180)

# This works for cst+ dst with 2 modes
sf=cst(rundir=[rundir, inv], damp=damp, model=None, 
       it_no=[enditer],lon_0=-180, cmap="su",
       map=True, plot=True, legend_show=False)#, savefigure=True)#, smax=2)
sf[0].stats.name = "Our inversion"

##### mzero
#cst_zero = join(rundir, invs, 'mzero.dat')
#sf += loadmodel(setup=setup, ifile=cst_zero, name='initial', damp='initial')       

##### IC + Mantle
net = '/net/home/talavera/eejit/splitting'   
ic =  'begtramp' #'GLW' #'tromp' "begtramp"

#main_mode = '01s00'
#main_mode = '02s00-07s02'
#main_mode = '03s00'
main_mode = '04s00'
#main_mode = '05s00'
#main_mode = '11s00-27s02'

#plot for RM paper
#cst_IC = '%s/%s/synthetics/IC/mcst-s20rts+crust+%s.dat'%(net,main_mode,ic)

#plot for IC paper
mode = '10s02-11s02'
#mode = '08s02-09s02'
cst_IC = '%s/%s/synthetics/IC/%s/mcst-s20rts+crust+%s.dat'%(net,main_mode,mode,ic) 
#cst_IC = '%s/%s/synthetics/IC/%s/mcst-%s.dat'%(net,main_mode,mode,ic) #ic only
#cst_IC = '%s/%s/synthetics/IC/mzero.dat'%(net,main_mode)

sf += loadmodel(setup=setup, ifile=cst_IC, 
#                name='%s'%ic, 
                name='Inner core model', 
                damp='0')       
 
model     = ['S20RTS'] 
#model     = ['AD', 'S20RTS']
#model     = ['S20RTS', 'AD']
#model = None
if model is not None:
    for x in model:
        sf += loadmodel(format=x, setup=setup, 
                        name ="S20RTS+CRUST5.1",   
#                        name=x, 
                        damp=0)  

p = sf.plot_map(
            kind_in_title=False, 
#            fig_size=[7,1.5],
            show_title=False,
            legend_show=False, 
            fig_abc='a',
            lon_0=-180,
            cmap="su",
#            vmax=7.2,
#            vmin=-5.4,
#            savefigure=True,
            )