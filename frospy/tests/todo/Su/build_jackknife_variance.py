#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 16:31:54 2018

@author: talavera
"""

import numpy as np
import os, fnmatch
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from frospy.util.read import read_modes
from frospy.preprocessing.spectrum import get_mode
from IPython.display import display, Math

if len(sys.argv) < 3:
    print "USAGE python bin_to_ascii.py damping mode(i.e 1S0)"
    sys.exit()


cwd = os.getcwd() # working directory
damp = float(sys.argv[1])
mode_name = sys.argv[2]

#damp = ('{0:.2e}').format(float('1000')) #damp='1.00e+03'
#mode_name='1S0'

jk_num  = len(list(filter(os.path.isdir, os.listdir(cwd))))-1
it_max  = len(fnmatch.filter(os.listdir('%s/1' %cwd), 'Iter_*'))
cst_num = int(open('%s/%d/Iter_1/cst_solver/new_model.damp_%.2e.cst' %(cwd, it_max, damp)).readline().rstrip())
cst_jk = []

# gettting cst coefficents from data
for i in range(1,jk_num+1):
    cst_file='%s/%d/Iter_%d/cst_solver/new_model.damp_%.2e.cst' %(cwd, i, it_max, damp)
    cst_jk.append(np.loadtxt(cst_file, skiprows=1))

# getting results using all events
cst_file= '%s/allevents/d%.0f/Iter_%d/cst_solver/new_model.damp_%.2e.cst' % (cwd, damp, it_max, damp)
cst_all = np.loadtxt(cst_file, skiprows=1)

# calculating variance
cst_var = [] # var = mean(abs(x - x.mean())**2).
for i in range(0,cst_num):
    cst_var.append(np.var(np.array(cst_jk).take([i], axis=1), ddof=cst_num))
    print r'cst%d = %.4f +/- %.1e' %(i, cst_all[i], cst_var[i])
    #display(Math(r'c_{st}^{%d} = %.4f \pm %.6f' %(i, cst_all[i], cst_var[i])))

# getting data from PREM
modedata = read_modes()
mode = get_mode([mode_name], modedata) #mode_sorted = sorted(mode.items(), key=lambda x: x[1]['freq'])
Q0 = mode[mode_name]['Q']; f0 = mode[mode_name]['freq']*1e3
fc = []; fc_var = []; Q = []; Q_var = []

fc.append( f0 + (4*np.pi)**-0.5 * cst_all[0] )
fc_var.append( fc[0] - (f0 + (4*np.pi)**-0.5 * cst_var[0]) )
Q.append( 0.5*fc[0] / (0.5*f0/Q0 + (4*np.pi)**-0.5 * cst_all[1]) )
Q_var.append( Q[0] - (0.5*fc[0] / (0.5*f0/Q0 + (4*np.pi)**-0.5 * cst_var[1])) )

print ''
print r'fc = %.4f +/- %.4f' %(fc[0], fc_var[0])
print r' Q = %.4f +/- %.4f' %(Q[0], Q_var[0])
#display(Math(r'f_{c} = %.4f \pm %.4f' %(fc[0], fc_var[0])))
#display(Math(r'Q = %.4f \pm %.4f' %(Q[0], Q_var[0])))

#plotting variance 
#x = np.linspace(1,cst_num,cst_num)
#ax = plt.subplot(1, 1, 1)
#ax.errorbar(x, cst_all, cst_var, marker='^', c='r') # linestyle='None', 
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.set_ylabel(r'$c_{st}$ ($\mu$Hz)')
#ax.set_xlabel('index')
#mf = ticker.ScalarFormatter(useMathText=True)
#mf.set_powerlimits((-2,2))
#plt.gca().yaxis.set_major_formatter(mf)
#plt.show()

