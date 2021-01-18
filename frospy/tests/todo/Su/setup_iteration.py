
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 09:55:09 2018

@author: talavera
"""
#--------------------------------------------------------------
# 1 mode inversion setup

import glob, os
import numpy as np
from os.path import basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.core.setup.settings import Setup
from frospy.core.setup.builder import build_inversion, build_array_inversion
from frospy.core.setup.create import rundir,dorun_addATA, segments, arraydirs
#from frospy.setup.create import dorun_pbs_header

main_mode  = '00s00'
event ='%sevents'%'Qc'
cpl   ='c00+d00' #c00+d00 cst+d00
mzero ='c00=1' #romanowicz cst+d00_best_result
damp = [0.0]
damp.extend(np.logspace(-4, -1, 4).tolist())
#damp=np.logspace(-4, 0, 5).tolist()
#damp=np.logspace(-4, -1, 4).tolist()

mdir = '/scratch/talavera/splitting'
#sdir = '/quanta1/home/talavera/splitting/%s/segments'%main_mode
sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-Qc'% (main_mode,
                                                                 main_mode)
#sdir2 = '/quanta1/home/talavera/splitting/%s/segments/13s02-all' % main_mode

input = OrderedDict() #AttribDict()
input['pbsjob'] = main_mode
input['rundir'] = os.path.join(mdir, main_mode, event, cpl, mzero)
input['bindir'] = '/quanta1/home/talavera/bin'
input['nmpydir'] = '/quanta1/home/talavera/codes/nmPy'
input['startmodel'] = 'PREM'
input['segmentsdir'] = sdir

#input['segmentsdir'] = AttribDict({'01s00': sdir, '04s02-00s10': sdir2})
input['damping'] = damp
input['iterations'] = [1, 10]
input['modes'] = OrderedDict([("00s00", 0),("02s02", 0)])#, ("11s02", 4)])#,
#                              ("10s16", 0), ("11s13", 0)])
#input['modes_cc'] = OrderedDict([("06s00-16s02",[0, 0])])
setup=Setup('CST', input)

#setup = build_inversion(input, 'CST', remove_existing=True)
setup = build_array_inversion(input, 'CST', remove_existing=True)

#--------------------------------------------------------------
# 2 mode inversion setup
import glob, os
import numpy as np
from os.path import basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.core.setup.settings import Setup
from frospy.core.setup.builder import build_inversion, build_array_inversion
from frospy.core.setup.create import rundir,dorun_addATA, segments, arraydirs

main_mode  = '10s02-11s02'
event ='%sevents'%'all'
cpl   ='cstC+dst' #c00+d00 cst+d00
mzero ='test' #romanowicz cst+d00_best_result
#damp = [0.0]
#damp.extend(np.logspace(-4, 0, 5).tolist())
damp=np.logspace(-4, 0, 5).tolist()
#damp=np.logspace(-4, -1, 4).tolist()

mdir = '/scratch/talavera/splitting'
sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-all'% (main_mode,
                                                                 main_mode)
#sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-good-AD'% (main_mode,
#                                                                 '02s00')
#sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-deep' % (main_mode,
#                                                                 '05s00')
#sdir2 = '/quanta1/home/talavera/splitting/%s/segments/13s02-all' % main_mode

input = OrderedDict() #AttribDict()
input['pbsjob'] = main_mode
input['rundir'] = os.path.join(mdir, main_mode, event, cpl, mzero)
input['bindir'] = '/quanta1/home/talavera/bin'
input['nmpydir'] = '/quanta1/home/talavera/codes/nmPy'
input['startmodel'] = 'PREM'
input['segmentsdir'] = sdir
#input['segmentsdir'] = AttribDict({'05s00': sdir, '13s02': sdir2})
input['damping'] = damp
input['iterations'] = [1, 10]
input['modes'] = OrderedDict([("10s02", 4), ("11s02", 4)])
input['modes_cc'] = OrderedDict([("10s02-11s02",[0, 4])])
input['modes_sc_dst'] = OrderedDict([("10s02", 2), ("11s02", 2)])
setup=Setup('CST', input)

#setup = build_inversion(input, 'CST', remove_existing=True)
setup = build_array_inversion(input, 'CST', remove_existing=True)

#--------------------------------------------------------------
# 3 mode inversion setup
import glob, os
import numpy as np
from os.path import basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.core.setup.settings import Setup
from frospy.core.setup.builder import build_inversion, build_array_inversion
from frospy.core.setup.create import rundir,dorun_addATA, segments, arraydirs

main_mode  = '04s00'
othr_mode = "10s02-11s02" 
event ='%sevents'%'all'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='good' #romanowicz cst+d00_best_result
#damp = [0.001]
#damp.extend(np.logspace(-4, 0, 5).tolist())
damp=np.logspace(-4, 0, 5).tolist()
#damp = [0.0001, 0.001, 0.01, 0.05, 0.1, 1.0]

mdir = '/scratch/talavera/splitting'
sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-good' % (main_mode,
                                                                    main_mode)
sdir2 = '/quanta1/home/talavera/splitting/%s/segments/%s-all' % (main_mode,
                                                                 othr_mode)

allmodes = "%s-%s" % (main_mode, othr_mode)

input = OrderedDict() #AttribDict()
input['pbsjob'] = allmodes
input['rundir'] = os.path.join(mdir, allmodes, event, cpl, mzero)
input['bindir'] = '/quanta1/home/talavera/bin'
input['nmpydir'] = '/quanta1/home/talavera/codes/nmPy'
input['startmodel'] = 'PREM'
#input['segmentsdir'] = sdir
input['segmentsdir'] = AttribDict({main_mode: sdir, othr_mode: sdir2})
input['damping'] = damp
input['iterations'] = [1, 10]
input['modes'] = OrderedDict([("04s00", 0), ("10s02", 4), 
                              ("11s02", 4)])#, ("00t11", 0)])
input['modes_cc'] = OrderedDict([("04s00-10s02",[2, 2]), 
                                 ("04s00-11s02",[2, 2]),
                                 ("10s02-11s02",[0, 4])])
#setup=Setup('CST', input)

#setup = build_inversion(input, 'CST', remove_existing=True)
setup = build_array_inversion(input, 'CST', remove_existing=True)

# -------------------------------------------------------------------
# dst inversion
import glob, os
import numpy as np
from os.path import basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.core.setup.settings import Setup, modes_sanity_check
from frospy.core.setup.builder import build_inversion, build_array_inversion
from frospy.core.setup.create import rundir,dorun_addATA, segments, arraydirs
#from frospy.setup.create import dorun_pbs_header

main_mode="10s02-11s02"
n = "10"; l='02'
deg_c=int(l)*2
if deg_c > 20:    
    deg_c=20
else:
    deg_c=deg_c
if deg_c > 6:    
    deg_d=6
else:
    deg_d=deg_c
    
R=0.2
modein = '%sS%s'%(int(n),int(l))
mode  = '%ss%s'%(n,l)
#modein = '%sT%s'%(int(n),int(l))
#mode  = '%st%s'%(n,l)
#cpl   ='cst' #c00+d00 cst+d00
cpl   ='cst+dst' #c00+d00 cst+d00
mzero ='c20' #romanowicz cst+d00_best_result
#mzero ='R=%s_TR'%R #romanowicz cst+d00_best_result
#degree = 'cst%s'%(deg_c)
degree = 'cst%s+dst%s'%(deg_c,deg_d)
#damp = [0.001]
damp=np.logspace(-4, -1, 4).tolist()

mdir = '/scratch/talavera/splitting'
#sdir = '/quanta1/home/talavera/splitting/%s/segments/good_80h' % (mode)
#sdir = '/quanta1/home/talavera/splitting/%s/segments/AD' % (mode)
#sdir2 = '/quanta1/home/simons/splitting/modes/%s/segments_T' % (mode)
#sdir = '/quanta1/home/talavera/splitting/%s/segments/all/segments_Z' % (mode)
#sdir2 = '/quanta1/home/talavera/splitting/%s/segments/all/segments_T'% (mode)
#sdir3 = '/quanta1/home/talavera/splitting/%s/segments/all/segments_R'% (mode)
sdir = '/quanta1/home/talavera/splitting/%s/segments/%s-all'% (main_mode,
                                                                 main_mode)

input = OrderedDict() #AttribDict()
input['pbsjob'] = mode
input['rundir'] = os.path.join(mdir, mode, cpl, mzero, degree)
input['bindir'] = '/quanta1/home/talavera/bin'
input['nmpydir'] = '/quanta1/home/talavera/codes/nmPy'
input['startmodel'] = 'PREM'
#input['startmodel'] = 'S20RTS'
#input['dst_R'] = R    
#input['startmodel'] = '%s/05s03/cst/prem/cst6/inversion_out/mnew-it10-d1.00e-04.dat'%mdir
#input['segmentsdir'] = sdir
#input['segmentsdir'] = AttribDict({'VHZ': sdir, 'VHT': sdir2, 'VHR': sdir3})
#input['segmentsdir'] = AttribDict({'VHT': sdir2, 'VHR': sdir3})
input['segmentsdir'] = AttribDict({'VHZ': sdir})
input['damping'] = damp
input['iterations'] = [1, 10]
input['modes'] = OrderedDict([("10s02", 4), ("11s02", 0)])
input['modes_sc_dst'] = OrderedDict([("10s02", 4), ("11s02", 0)])
# input['modes'] = OrderedDict([(modein, deg_c), ("00t%s"%(int(l)+1), 0)])
#input['modes'] = OrderedDict([(modein, deg_c)])
#input['modes_sc_dst'] = OrderedDict([(modein, deg_d)])
#input['modes_sc_dst'] = OrderedDict([(modein, deg_d), ("00t%s"%(int(l)+1), 0)])
setup=Setup('CST', input)

#setup = build_inversion(input, 'CST', remove_existing=True)
setup = build_array_inversion(input, 'CST', remove_existing=True)
#-------------------------------------------------------------------
# Uncertainties
from frospy.core.setup.builder import build_uncertainties
import frospy.core.setup.create as create
from frospy.core.setup.settings import read as read_setup
from frospy.core.setup.builder import build_array_inversion


#mdir = '/quanta1/home/talavera/splitting'
mdir = '/scratch/talavera/splitting'
#mode  = '00s00'
#event ='%sevents'%'Qc'
#cpl   ='c00+d00' #c00+d00 cst+d00
#mzero ='c00=1' #romanowicz cst+d00_best_result

n = "00"; l='06'
deg_c=int(l)*2
if deg_c > 20:    
    deg_c=26
else:
    deg_c=deg_c
if deg_c > 6:    
    deg_d=6
else:
    deg_d=deg_c
#modein = '%sS%s'%(int(n),int(l))
#mode  = '%ss%s'%(n,l)
modein = '%sT%s'%(int(n),int(l))
mode  = '%st%s'%(n,l)
cpl   ='cst+dst' #c00+d00 cst+d00
mzero ='prem_TR' #romanowicz cst+d00_best_result
degree = 'cst%s+dst%s'%(deg_c,deg_d)
damp = 1.00e-02

#run_dir = os.path.join(mdir, mode, event, cpl, mzero)
run_dir = os.path.join(mdir, mode, cpl, mzero, degree)
pickle_dir= os.path.join(run_dir, 'setup.json')
setup     = read_setup(pickle_dir)
setup.walltime='360:00:00'
create.rundir(setup, remove_existing=False, update_only=True)
setup = build_array_inversion(setup, 'CST', remove_existing=False,
                              scripts_only=True)
setup = build_uncertainties(run_dir, damp, 
                            remove_existing=False, 
                            boot_strap=False, 
                            jack_knive=True,
                            allevents_subset_ratio=0.7,
#                            N_of_subsets=20,
                            )
print(run_dir)
from frospy.core.setup.settings import read as read_setup
from frospy.postprocessing.uncertainties import uncertainties_calculation
from os.path import join

#cwd   = '/scratch/talavera/splitting'
#cwd   = '/Users/sujania/Documents/PhD/UU.eejit/splitting'
cwd = '/quanta1/home/talavera/splitting'
#cwd = '/net/home/talavera/eejit/splitting'
mode  = '10s02-11s02'
event ='%sevents'%'all'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='c00' #romanowicz cst+d00_best_result
inv=1 
damp = 0.01

pickle_dir= join(cwd,mode,event,cpl,mzero, 'setup.json')
setup     = read_setup(pickle_dir)
db_path = join(cwd,mode,event,cpl,mzero,'uncertainties','%s.sqlite3'%mode)

SF = uncertainties_calculation(setup, damp, 
                               db=True, 
                               db_path=db_path, 
                               model='data',
                               verbose=True,
                               jack_knive=True,)


#-------------------------------------------------------------------
# Update old innversion
from frospy.core.setup.settings import read, update_setup
from frospy.core.setup.builder import build_array_inversion
from frospy.core.setup import create
from os.path import join

mdir = '/quanta1/home/talavera/splitting'
mode  = '10s02-11s02'
event ='%sevents'%'all'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='c00' #romanowicz cst+d00_best_result

run_dir = os.path.join(mdir, mode, event, cpl, mzero)

setup = read(join(run_dir, 'setup.json'))
setup = update_setup(setup, 'quanta')
create.rundir(setup, remove_existing=False, update_only=True)
setup = build_array_inversion(setup, 'CST', remove_existing=False,
                              scripts_only=True)