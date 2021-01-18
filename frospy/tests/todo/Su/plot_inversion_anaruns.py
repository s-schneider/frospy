## damping
import matplotlib.pyplot as plt
from frospy.postprocessing.AnalyseRuns.plot import misfits_cst, listdir, cst_map
from frospy.core.splittingfunc import Set
from frospy.core.splittingfunc import loadmodel # python2 branch
#from frospy.splitting.load import loadmodel # python3 master branch
from frospy.postprocessing.misfit import plot_Qf
from frospy.util.read import read_modes_in
import numpy as np
import os

cwd   = os.getcwd()
mode  = os.path.basename(cwd)
event ='%sevents'%'all'
cpl   ='cst+d00' #cst+d00
mzero ='cst_01s00' #romanowicz cst+d00_best_result
damp  = 'd0.00001'#d0.01
#damp  = 'best_result'#d0.01
cst_deg = range(14,18+2,2)

#run_dir   = '%s/%s/%s/%s'%(cwd,event,cpl,mzero)
run_dir = '//nfs/stig/talavera/radial-inversion/05s00/deepevents/cst+dst/cst+d00_best_result'
title     = '%s with %s events in %s-cpl, mzero=%s'%(mode, event,cpl,mzero)
dampings  = [10000, 10000, 1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
#dampings  = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
startiter = 1
enditer   = 10
iters     = range(startiter, enditer+1)

if damp != '':
# for a specific damping
    run_dir='%s/%s'%(run_dir,damp)
    misfits_cst(run_dir=run_dir, comp="Z", dampings=dampings, iterations=iters,
                title='%s, with d=%s'%(title,damp), cst_deg=cst_deg, 
                plot_coeff=True, models=['S20RTS'])
    cst_map(run_dir=run_dir, iteration=enditer, 
            models=None, smax=max(cst_deg)) 
#            title='%s, with d=%s'%(title,damp))
#    cst_map(run_dir=run_dir, iteration=enditer, models=['S20RTS'], smax=2) 
else:
    runs = listdir(run_dir, 'd'); runs.sort()
    misfits_cst(run_dir=run_dir, comp="Z", dampings=dampings, 
                iterations=iters, title=title, cst_deg=cst_deg, 
                plot_coeff=True, models=['S20RTS'])
#plt.close('all')

##----- tests
#run_dir='%s/%s'%(run_dir,damp)
##run_dir=os.getcwd()
#cst_file = '%s/Iter_%s/cst_solver/mcst.dat' % (run_dir, enditer)
##cst_file = '%s/config/mcst_zero.dat' % (run_dir)
#modes_dir = '%s/config' % run_dir
##
#splf = loadmodel(cst_file, modesin_dir=modes_dir)
#splf.plot_map(vmin=-7,vmax=7)
#

#from frospy.core.splittingfunc.read import  read_cst_S20RTS
#modesin, modes_ccin = read_modes_in(modes_dir)
#cst=read_cst_S20RTS(modesin, modes_ccin)

#splf = loadmodel(cst_file, modesin_dir=modes_dir)
#splf.plot_map(smax=12)
#splf = loadmodel(format='S20RTS', modesin_dir=modes_dir)
#splf.plot_map(smax=12, filename='S20RTS')

#splf = loadmodel(cst_file, modesin_dir=modes_dir)
#splf.plot(smax=4)
#splf = loadmodel(format='S20RTS', modesin_dir=modes_dir)
#splf.plot(smax=12, filename='S20RTS')
#from frospy.core.splittingfunc.read import read_cst
#cst=read_cst(cst_file, modes_dir=modes_dir)
    
from frospy.plot.nmplt import plot_Mmatrix
run_dir = '//nfs/stig/talavera/mantle-inversion/04s02/04s02-00s10-00t11'
enditer   = 9
damp = '0.001'
modes_dir = '%s/d%s/config' % (run_dir,damp)
matrix_dir = '%s/d%s/Iter_%s/mdcpl/matrix.dat'%(run_dir, damp, enditer)
omega_dir = '%s/d%s/Iter_%s/matdiag/omega.dat'%(run_dir, damp, enditer)
plot_Mmatrix(matrix_dir, omega_dir, modes_dir)