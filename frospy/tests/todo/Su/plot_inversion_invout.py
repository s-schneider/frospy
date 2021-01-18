#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 11:42:56 2018

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
from frospy.core.setup.settings import update_setup_pickle
from frospy.util.read import read_pickle
#---------------------------------------------------------------------------
# 1 modes
cwd   = '/net/home/talavera/eejit/splitting'
#cwd   = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
mode  = '00s06'
event ='cst+dst'#'%sevents'%'all'
cpl   ='syn_R=0_TRZ' #R=0.2_AD prem_TRZ
mzero ='cst12+dst12' 
#event ='cst'#'%sevents'%'all'#
#cpl   ='syn_R=-0.2' #c00+d00 cst+d00 prem_TRZ
#mzero ='cst20' 
#event ='%sevents'%'Qc'#'cst'#
#cpl   ='c00+d00' #c00+d00 cst+d00
#mzero ='c00=1' #romanowicz cst+d00_best_result
inv=1
lcut=6
R=-0.2

#pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.pickle')
pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.json')
title     = '%s with %s events in %s, mzero=%s'%(mode, event,cpl,mzero)
setup     = read_setup(pickle_dir)
rundir    = join(cwd,mode,event,cpl,mzero)
dampings  = setup.damping
startiter = setup.iterations[0]
enditer   = 5#setup.iterations[1]
iters     = range(startiter, enditer+1)
damp      = 0.01
#model     = ['S20RTS', 'AD']
#model     = ['S20RTS']
model     = None

#["VHZ", "VHR", "VHT"]
inv_summary(schema=["VHZ", "VHR", "VHT"],rundir=[rundir, inv], title=title,
#inv_summary(schema=["VHZ", "VHT"],rundir=[rundir, inv], title=title,
#inv_summary(schema=["VHR", "VHT"],rundir=[rundir, inv], title=title,
#inv_summary(schema=["VHZ"],rundir=[rundir, inv], title=title,
#inv_summary(schema=['05s00','13s02'], rundir=[rundir,inv], title=title, 
#inv_summary(schema=["T", "R"],rundir=[rundir, inv], title=title,
#inv_summary(schema='Set',rundir=[rundir, inv], title=title, 
            print_summary=True, plot=True, lcurve=True,colormap='rainbow',)#,
#            damping=[1.00e-00,1.00e-01])

cst(rundir=[rundir, inv], model=model, it_no=[enditer], cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp="all",
    smax=lcut)#,damp=[1.00e-01, 1.00e-02,1.00e-03])
    #, verbose=True)
cst(rundir=[rundir, inv], model=model, it_no="all",smax=lcut, cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp=damp)#, verbose=True)
sf=cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=lcut,R=R)#, savefigure=True)

cst(rundir=[rundir, inv], model=model, it_no=[enditer], smax=2, cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp="all",)#, verbose=True) 

#cst(rundir=[rundir, inv], model=model, it_no=[enditer], cmap='rainbow',
#    map=False, plot=True, plot_mzero=True, damp="all",
#    smax=2)#,damp=[1.00e-01, 1.00e-02,1.00e-03])
#    #, verbose=True)s
#       
#cst(rundir=[rundir, inv], model=model, it_no="all",smax=2, cmap='rainbow',
#    map=False, plot=True, plot_mzero=True, damp=damp)#, verbose=True)
#cst(rundir=[rundir, inv], model=model, it_no="all",smax=8, cmap='rainbow',
#    map=False, plot=True, plot_mzero=True, damp=damp)#, verbose=True)


cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=2,R=R)#, savefigure=True)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=4,R=R)#, savefigure=True)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=6,R=R)#, savefigure=True)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=8,R=R)#, savefigure=True)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, legend_show=False, smax=10,R=R,)#cmap='su')#, savefigure=True)
#
#setup.iterations[1] = 5
#setup.modes_sc = OrderedDict([("0S6", 12)])
#setup.model_size = 182
#setup = update_setup_pickle(setup)
#setup.write(join(cwd,mode,event,cpl,mzero))#, overwrite=True)
#---------------------------------------------------------------------------
# 2 modes

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
from frospy.core.setup.settings import update_setup_pickle
from frospy.util.read import read_pickle

cwd   = '/net/home/talavera/eejit/splitting'
#cwd   = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
#mode  = '05s00-13s02'
#mode  = '08s02-09s02'
#mode  = '10s02-11s02'
mode  = '00s00-00s05'
#event ='%sevents'%'deep'
event ='%sevents'%'Qc'
cpl   ='cst+d00+cc' #c00+d00 cst+d00 cst+d00+cc cstC+dst
#mzero ='cst_05s00_13s02' #romanowicz cst+d00_best_result
mzero ='c00=1'
inv=4# 4, 5, 6, 9, [11, (16)], [12, (13, 14, 15)]

pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.json')
setup     = read_setup(pickle_dir)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
rundir    = join(cwd,mode,event,cpl,mzero)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
dampings  = setup.damping
startiter = setup.iterations[0]
enditer   = 10#setup.iterations[1]
iters     = range(startiter, enditer+1)
damp      = 0.001
#model     = ['S20RTS', 'AD']
model     = ['S20RTS']
#model     = None

inv_summary(schema=["00s05", "00s00"],rundir=[rundir, inv], colormap='rainbow',
#inv_summary(schema='Set',rundir=[rundir, inv], colormap='rainbow',
            title=title, bins=30,lcurve=True,
            print_summary=True, plot=True, damping="all")#
#            damping=[1.00e-01, 1.00e-02, 1.00e-03])

cst(rundir=[rundir, inv], model=model, it_no=[enditer], cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp="all")#
#    damp=[1.00e+00,1.00e-01, 1.00e-02, 5.00e-02])#, verbose=True)

cst(rundir=[rundir, inv], model=model, it_no="all", cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp=damp)#, verbose=True)
sf=cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer],lon_0=-180, 
    map=True, plot=True, legend_show=False)#, savefigure=True)#, smax=2)
#setup.modes_sc_dst['10s02']=2
#setup.modes_sc_dst['11s02']=2
#setup.model_size=58
#setup = update_setup_pickle(setup)
#setup.write(join(cwd,mode,event,cpl,mzero))#, overwrite=True)
#### With IC and other models
# best result
sf = Set()
if inv == 1:
    invs = 'inversion_out'
else:
    invs = 'inversion_out_%s'%inv
    
cst_files = glob.glob('%s/%s/mnew-it%s-d%.2e.dat' % 
                      (rundir, invs, enditer, damp))
for i, f in enumerate(cst_files):
    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
    sf += loadmodel(setup=setup, ifile = f, 
#                    name ="Our measurements", 
                    name ="it=%s, d=%s"%(str(it),damp), 
                    damp = damp)

# mzero
cst_zero = join(rundir, invs, 'mzero.dat')
#sf += loadmodel(setup=setup, ifile=cst_zero, name='initial', damp='initial')       

# IC + Mantle
ic =  's20rts+crust+tromp' #'GLW' #
net = '/net/home/talavera/radial-inversion'
#cst_IC = '%s/04s00/synthetics/IC/%s/mcst-%s.dat'%(net,mode,ic)   
cst_IC = '%s/11s00/synthetics/IC/mcst-%s.dat'%(net,ic) 
sf += loadmodel(setup=setup, ifile=cst_IC, 
                name='%s'%ic, 
#                name='Inner core model', 
                damp='0')       
 
model     = ['S20RTS'] 
#model     = ['AD', 'S20RTS']
#model     = ['S20RTS', 'AD']
if model is not None:
    for x in model:
        sf += loadmodel(format=x, setup=setup, 
                        name ="S20RTS+CRUST5.1",   
#                        name=x, 
                        damp=0)  

sf.plot_map()#kind_in_title=False)
#---------------------------------------------------------------------------
# CC cpl
cwd   = '/net/home/talavera/eejit/splitting'
#cwd   = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
mode  = '03s00-08s02-09s02'#
#mode  = '04s00-10s02-11s02'#
#mode  = '01s00-04s02-00s10'#
event ='%sevents'%'good'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='08s02-09s02'#'08s02-09s02' #romanowicz cst+d00_best_result
#mzero ='10s02-11s02'
inv   = 1

pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.json')
setup     = read_setup(pickle_dir)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
rundir    = join(cwd,mode,event,cpl,mzero)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
dampings  = setup.damping
startiter = setup.iterations[0]
enditer   = 15#setup.iterations[1]
iters     = range(startiter, enditer+1)
damp      = 0.05
model     = ['S20RTS']#, 'AD']

#inv_summary(schema=['04s00','10s02-11s02'], fig_size=(10, 10),
#inv_summary(schema=['01s00','04s02-00s10'], fig_size=(10, 10),
inv_summary(schema=['03s00','08s02-09s02'], fig_size=(10, 10),
            rundir=[rundir, inv], title=title, #bins=30,
            print_summary=True, plot=True,lcurve=True,)
#            damping=[1.00e-00,1.00e-01, 1.00e-02, 5.00e-02])#, 1.00e-04])

sf=cst(rundir=[rundir, inv], model=model, it_no=[enditer], cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp="all")#
#    damp=[1.00e-00,1.00e-01, 1.00e-02, 5.00e-02])#, verbose=True)
#            damp=[1.00e-00,1.00e-01, 1.00e-02, 5.00e-03, 1.00e-03])#, 1.00e-04])

cst(rundir=[rundir, inv], model=model, it_no="all", cmap='rainbow',
    map=False, plot=True, plot_mzero=True, damp=damp)#, smax=4)#, verbose=True)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True)#, smax=2)

#### With IC and other models
# best result
sf = Set()
cst_files = glob.glob('%s/%s/mnew-it%s-d%.2e.dat' % 
                      (rundir, 'inversion_out', enditer, damp))
for i, f in enumerate(cst_files):
    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
    sf += loadmodel(setup=setup, ifile = f, 
                    name ="it=%s, d=%s"%(str(it),damp), damp = damp)

## mzero
#cst_zero = join(rundir, 'inversion_out', 'mzero.dat')
#sf += loadmodel(setup=setup, ifile=cst_zero, name='initial', damp='initial')       

# IC + Mantle
ic = 's20rts+crust+begtramp' #'GLW' # 
net = '/net/home/talavera/radial-inversion'
cst_IC = '%s/03s00/synthetics/IC/mcst-%s.dat'%(net,ic)   
#sf += loadmodel(setup=setup, ifile=cst_IC, name='%s'%ic, damp='0')            
  
#model     = None 
#model     = ['S20RTS'] 
#model     = ['AD', 'S20RTS']
model     = ['S20RTS', 'AD']
if model is not None:
    for x in model:
        sf += loadmodel(format=x, setup=setup, name=x, damp=0)  
        
sf.plot_map()

#---------------------------------------------------------------------------
# dst inversion
cwd    = '/net/home/talavera/eejit/splitting'
#cwd    = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
mode   = '01s04'
cpl    ='cst+dst' #c00+d00 cst+d00
mzero  ='prem_TRZ' #romanowicz cst+d00_best_result
degree = 'cst8+dst8'#'cst14+dst8'
inv    = 1

pickle_dir= join(cwd,mode,cpl,mzero,degree,'setup.json')
setup     = read_setup(pickle_dir)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode,cpl,mzero,degree)
rundir    = join(cwd,mode,cpl,mzero,degree)
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode,cpl,mzero,degree)
dampings  = setup.damping
startiter = setup.iterations[0]
enditer   = setup.iterations[1]
iters     = range(startiter, enditer+1)
damp      = 0.1
model     = ['S20RTS']#, 'AD']
R         = -0.2
lcut      = 2
#inv_summary(schema='Set',rundir=[rundir, 1], title=title, bins=30,
#            print_summary=True, plot=True)#,
##            damping=[1.00e-00, 1.00e-01, 1.00e-02])

cst(rundir=[rundir, inv], damp='all', model=model, it_no=[enditer], cmap='rainbow',
    map=False, plot=True, plot_mzero=True, smax=lcut, R=R)
cst(rundir=[rundir, inv], damp=damp, model=model, it_no=[enditer], 
    map=True, plot=True, R=R)#, smax=2)

#---------------------------------------------------------------------------
# All dampings
damp = 0.01
R    = -0.2
cst_files = glob.glob('%s/%s/mnew-it%s-d%.2e.dat' % 
                      (rundir, 'inversion_out', enditer, damp))
#cst_files = glob.glob('%s/%s/mnew-it%s-d*.dat' % 
#                      (rundir, 'inversion_out', enditer))
sort_human(cst_files)
cst_zero = join(rundir, 'mzero.dat')
sf = Set()
#sf += loadmodel(setup=setup, ifile=cst_zero, name='mzero', damp='initial')       
models = ['S20RTS']#, 'AD']

mode = "2S6" # correcting dst ellipticity
for i, f in enumerate(cst_files):
    damp = float(splitext(f)[-2].split('-d')[-1])    
    sf += loadmodel(setup=setup, ifile = f, name = "d=%s, it=%s" %(damp, 10), 
                    damp = damp)
    # correcting dst ellipticity
#    sf[i].dst[mode]["2"][0] = sf[i].dst[mode]["2"][0]/10
#    sf[i+1].dst[mode]["2"][0] = sf[i+1].dst[mode]["2"][0]/10
if models is not None:
    for x in models:
        if x is "S20RTS": 
            sf += loadmodel(format=x, setup=setup, name="%s, R=%s" % (x, R), 
                            damp=0, R=R)  
        else:
            sf += loadmodel(format=x, setup=setup, name=x, damp=0)  
#sf.plot()
sf.plot_map(smax=4, kind="cst")
sf.plot_map(smax=4, kind="dst")#, )
#sf[4].plot_map(savefig=True)
#sf[6].plot_map()
#plt.close('all')

#---------------------------------------------------------------------------
# All iterations
damp      = 0.1
cst_files = glob.glob('%s/%s/mnew-it*-d%.2e.dat' % 
                      (rundir, 'inversion_out', damp))
#cst_files = glob.glob('%s/%s/mnew-it%s-d*.dat' % 
#                      (rundir, 'inversion_out', enditer))
cst_zero = join(rundir, 'inversion_out', 'mzero.dat')
sort_human(cst_files)
models = ['S20RTS', 'AD']
sf = Set()
sf += loadmodel(setup=setup, ifile=cst_zero, name='initial', damp='initial')       
for i, f in enumerate(cst_files):
    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
    sf += loadmodel(setup=setup, ifile = f, name = str(it), damp = damp)
     
if models is not None:
    for x in models:
        sf += loadmodel(format=x, setup=setup, name=x, damp=0)  
sf.plot()
#sf[4].plot_map(savefig=True)
#sf[6].plot_map()
#plt.close('all')

#--------------------------------------------------------------------------
# Measurement Comparisson 2 modes
cwd   = '/net/home/talavera/eejit/splitting'
mode  = '01s00-04s02-00s10'
event ='%sevents'%'all'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='prem' #romanowicz cst+d00_best_result
rundir    = join(cwd,mode,event,cpl,mzero)
pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.pickle')
setup     = read_setup(pickle_dir)

c    = 0 # s20rts+crust+ic selection
ic   = 'tromp'
damp = 0.01

sf = Set()
cst_files = glob.glob('%s/%s/mnew-it10-d%.2e.dat' % 
                      (rundir, 'inversion_out', damp))
sort_human(cst_files)    
        
for i, f in enumerate(cst_files):              
    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
    sf += loadmodel(setup=setup, ifile = f, name = str(it), damp = damp)
 
models = None#['AD']#'S20RTS']     
if models is not None:
    for x in models:
        sf += loadmodel(format=x, setup=setup, name=x, damp=0)
        
radial = ['02s00','05s00','11s00']
innerc = ['07s02','13s02','27s02']
cdir1 = '/net/home/talavera/radial-inversion'
cdir2 = 'synthetics/IC/mcst-s20rts+crust+%s.dat'%ic
cst_files = [join(cdir1,m,cdir2) for m in radial]

for f in [cst_files[c]]:
    input2 = OrderedDict() #AttribDict()
    input2['rundir'] = cst_files[c]
    input2['modes'] = OrderedDict([(radial[c], 0), (innerc[c], 4)])
    input2['modes_cc'] = OrderedDict([('%s-%s'%(radial[c], innerc[c]),[2, 2])])
    input2['startmodel'] = 'PREM'
    setup2=Setup('CST', input2)
    sf += loadmodel(setup=setup2, ifile = f, 
                    name = 'S20RTS+CRUST+IC', damp = 'S20RTS+CRUST+IC')    

models = ['S20RTS']     
if models is not None:
    for x in models:
        sf += loadmodel(format=x, setup=setup, name=x, damp=0)
     
sf.plot_map()

#--------------------------------------------------------------------------
# Measurement Comparisson 3 modes
import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.core.setup.settings import Setup, get_mode_list
from frospy.core.splittingfunc import loadmodel, Set
from frospy.core.setup.settings import read as read_setup

sf = Set()
models = None#['AD']#'S20RTS']  
if models is not None:
    for x in models:
        sf += loadmodel(format=x, setup=setup, name=x, damp=0)
        
#ic = ['begtramp']
ic = ['begtramp', 'romanowicz', 'tromp', 'woodhouse']
#radial = ['03s00','04s00']
#innerc1 = ['08s02','10s02']
#innerc2 = ['09s02','11s02']
radial = ['04s00']
innerc1 = ['10s02']
innerc2 = ['11s02']
cdir1 = '/net/home/talavera/radial-inversion'
cdir2 = [['synthetics/IC/%s-%s/mcst-s20rts+crust+%s.dat'%(ic1,ic2,i) 
        for ic1, ic2 in zip(innerc1, innerc2)] for i in ic]
cst_files = [join(cdir1,m,c[i]) for i,m in enumerate(radial) for c in cdir2]

#i = 0
#for f in [cst_files[i]]:
for i,r in enumerate(radial):
    for f in cst_files: 
#        sf = Set()
        if f.find(r) > 0:
            m = f.split('+')[-1].split('.')[0]
            if m == 'romanowicz':
                name = 'Durek & Romanowicz'#1999
            if m == 'woodhouse':
                name = 'Woodhouse et al.'#1986
            if m == 'tromp':
                name = 'Tromp'#1993
            if m == 'begtramp':
                name = 'Beghein & Trampert'#2003
            
            input2 = OrderedDict()#AttribDict()
            input2['rundir'] = cst_files[i]
            input2['modes'] = OrderedDict([(innerc1[i], 4), 
                                          (innerc2[i], 4)])
            input2['modes_cc'] = OrderedDict(
                                    [('%s-%s'%(innerc1[i],innerc2[i]),[0, 4])])
            input2['startmodel'] = 'PREM'
            setup2=Setup('CST', input2)
            sf += loadmodel(setup=setup2, ifile = f, 
                            name = name, damp = 'S20RTS+CRUST+IC') 
            s = loadmodel(setup=setup2, ifile = f, 
                            name = name, damp = 'S20RTS+CRUST+IC') 

    #damp      = 0.01
    #cst_files = glob.glob('%s/%s/mnew-it10-d%.2e.dat' % 
    #                      (rundir, 'inversion_out', damp))
    #sort_human(cst_files)            
    #for i, f in enumerate(cst_files):
    #    it = int(splitext(f)[-2].split('-it')[-1].split('-')[0])
    #    sf += loadmodel(setup=setup, ifile = f, name = str(it), damp = damp)
        
        s.plot_map(kind_in_title=False, 
                #   fig_size=[7,1.5],
                   show_title=False,
                    legend_show=False, 
                    cmap='su',
                    lon_0=180.0,
                    weight="bold",
                    filename=m,
    #               vmax=50,
    #               vmin=-50,
                    savefig=True,)
    sf.plot_map(kind_in_title=False, 
            #   fig_size=[7,1.5],
#               show_title=False,
                legend_show=False, 
                cmap='su',
                lon_0=180.0,
                weight="normal",
#               vmax=50,
#               vmin=-50,
                savefigure=False,)
#---------------------------------------------------------------------------
# 2 modes cross-coupling
#radial = ['02s00','05s00','06s00','07s00','08s00','09s00','11s00']
#innerc = ['07s02','13s02','16s02','18s02','20s02','22s02','27s02']

ic = 'tromp'
radial = ['02s00','05s00','11s00']
innerc = ['07s02','13s02','27s02']
cdir1 = '/net/home/talavera/radial-inversion'
cdir2 = 'synthetics/IC/mcst-s20rts+crust+%s.dat'%ic
cst_files = [join(cdir1,m,cdir2) for m in radial]

sf = Set()
for i, f in enumerate(cst_files):
    input = OrderedDict() #AttribDict()
    input['rundir'] = cst_files[i]
    input['modes'] = OrderedDict([(radial[i], 0), (innerc[i], 4)])
    input['modes_cc'] = OrderedDict([('%s-%s'%(radial[i], innerc[i]),[2, 2])])
    input['startmodel'] = 'PREM'
    setup=Setup('CST', input)
    sf += loadmodel(setup=setup, ifile = f, 
                    name = 'S20RTS+CRUST+%s'%ic, damp = 'S20RTS+CRUST+IC')    
    sf[i].plot_map()


radial = ['03s00','04s00']
innerc1 = ['08s02','10s02']
innerc2 = ['09s02','11s02']
cdir1 = '/net/home/talavera/radial-inversion'
cdir2 = 'synthetics/IC/mcst-s20rts+crust+%s.dat'%ic
cst_files = [join(cdir1,m,cdir2) for m in radial]

sf = Set()
for i, f in enumerate(cst_files):
    input = OrderedDict()#AttribDict()
    input['rundir'] = cst_files[i]
#    input['modes'] = {radial[i]: 0, innerc1[i]:4, innerc2[i]:4}
    input['modes'] = OrderedDict([(radial[i], 0), 
                                  (innerc1[i], 4), 
                                  (innerc2[i], 4)])
    input['modes_cc'] = OrderedDict([('%s-%s'%(radial[i], innerc1[i]),[2, 2]),
                                     ('%s-%s'%(radial[i], innerc2[i]),[2, 2]),
                                     ('%s-%s'%(innerc1[i],innerc2[i]),[0, 4])])
    input['startmodel'] = 'PREM'
    setup=Setup('CST', input)
    sf += loadmodel(setup=setup, ifile = f, 
                    name = 'S20RTS+CRUST+%s'%ic, damp = 'S20RTS+CRUST+IC')    
    sf[i].plot_map()
#    sf[i].plot()
#---------------------------------------------------------------------------
from frospy.core.database.query import get_db_tables
get_db_tables(os.path.join(rundir, 'inversion_out/inversion.db'))

from frospy.core.database.query import db_query
damp = 0.0001
sset = "01s00" # "01s00"
misfit = db_query(os.path.join(rundir, 'inversion_out/inversion.db'), sset,
                               condition_dict={'iter': enditer,
                                               'damp': damp},
                               select='misfit')
misfit2 = misfit

sset = "04s02-00s10" # 
misfit = db_query(os.path.join(rundir, 'inversion_out/inversion.db'), sset,
                                condition_dict={'iter': enditer,
                                                'damp': damp},
                               select='misfit')
print(np.average(misfit2),np.average(misfit))
misfit2.extend(misfit)
print(np.average(misfit2))

from frospy.util.read import (
    read_modes_in, get_mode_names, get_mode_deg, read_pickle
                            )
modesin, modes_ccin = setup.get_modes_in()
cst, dst, modes_sc, modes_cc = read_cst(setup=setup, cfile='S20RTS')
sc_modes, cc_modes = get_mode_names(modesin, modes_ccin)
sf = loadmodel(format="S20RTS", setup=setup, name="S20RTS", damp=0)
sf.plot_map()
#---------------------------------------------------------------------------

from frospy.plot.nmplt import plot_Mmatrix
run_dir = '//nfs/stig/talavera/mantle-inversion/04s02/04s02-00s10-00t11'
enditer   = 9
damp = '0.001'
modes_dir = '%s/d%s/config' % (run_dir,damp)
matrix_dir = '%s/d%s/Iter_%s/mdcpl/matrix.dat'%(run_dir, damp, enditer)
omega_dir = '%s/d%s/Iter_%s/matdiag/omega.dat'%(run_dir, damp, enditer)
plot_Mmatrix(matrix_dir, omega_dir, modes_dir)