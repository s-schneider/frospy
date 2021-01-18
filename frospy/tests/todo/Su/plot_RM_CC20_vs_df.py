#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 18:34:34 2019

@author: talavera
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os

#for m in ['06s00', '07s00', '08s00', '09s00', "11s00"]:
#    cwd = '/net/home/talavera/radial-inversion/%s/synthetics/LM'%m
#    df = []
#    with open("%s/summary.out"%cwd, "r") as file:
#        for line in file:
#            df.append(line.split())
#  
#    df = zip(*df)    
#    c20 = []
#    dfreq = []
#    for c, f in zip(df[0][:-6],df[1][:-6]):
#        c20.append(int(c))
#        dfreq.append(float(f))
#    c20.append(0)
#    dfreq.append(float(df[1][-4]))
#    if m is not "11s00":
#        plt.axhline(y=float(df[1][-6]), color='r', linestyle='-', label='REM') #REM
#    plt.axhline(y=float(df[1][-5]), color='b', linestyle='--', label='begtramp') #begtramp
#    plt.axhline(y=float(df[1][-3]), color='c', linestyle='-', label='romanowicz') #romanowicz
#    plt.axhline(y=float(df[1][-2]), color='g', linestyle='--', label='tromp') #tromp
#    plt.axhline(y=float(df[1][-1]), color='m', linestyle='-', label='woodhouse') #woodhouse
#    plt.plot(c20, dfreq, 'o', label=r'$c_{20}^{CC} \ne 0 $')
#    
#    xmin=35; xmax=105
#    if m=='06s00': plt.xlim(-xmin,xmax); plt.title(r'%s, r+e: $c_{20}^{CC}=$%d'%(m,18))
#    if m=='07s00': plt.xlim(-xmax,xmin); plt.title(r'%s, r+e: $c_{20}^{CC}=$%d'%(m,-22))
#    if m=='08s00': plt.xlim(-xmin,xmax); plt.title(r'%s, r+e: $c_{20}^{CC}=$%d'%(m,25))
#    if m=='09s00': plt.xlim(-xmax,xmin); plt.title(r'%s, r+e: $c_{20}^{CC}=$%d'%(m,-29))
#    if m=='11s00': plt.xlim(-xmax,xmin); plt.title(r'%s, r+e: $c_{20}^{CC}=$%d'%(m,-35))
#    plt.legend()
#    plt.show()
    
model = ["s20rts+crust","ellip","begtramp","tromp","woodhouse"]#,"romanowicz"]
ellip = [6.71, -17.68, 17.97, -21.80, 24.69, -28.71, -34.54]
PREM= [[2.51048, 1242.24], [4.88417, 920.81], [5.74025, 913.242], 
       [6.58071,881.057], [7.42413,851.789], [8.26264,840.336], 
       [9.88792,831.947]]
REM = [[2.50795, 1250.00], [4.88840, 1129.94], [5.7422,1063.83],
       [6.5850,892.86],[7.4290,1020.41],[8.2692,1250.00]]

SC = [[2.50816, 1015.92], [4.88837, 1020.04], [5.7423,1067.11],[6.5844,723.26],
        [7.4302,1221.95],[8.2701,946.62],[9.8865,838.88]]

# 2s0, 5s0
CCcoef = [5.8557348, -8.6669769]
CC = [[2.508023, 1201.44], [4.888208, 1000.32]] 

c = 0
for mode in ['02s00', '05s00']:#, '06s00', '07s00', '08s00', '09s00', "11s00"]:
#    f, (ax1, ax2) = plt.subplots(1, 2)
    f, (ax1) = plt.subplots(1, 1)
    cwd = '/net/home/talavera/radial-inversion/%s/synthetics/LM'%mode
    omegaR = []
    for c20 in np.arange(-100,110, 10):
        ofile = np.array([c20])
        ofile = np.append(ofile,
                          np.loadtxt("%s/omega.CC20=%s.dat"%(cwd,c20))[0])
        omegaR.append(ofile) 
        
#        ic = np.loadtxt("%s/omega.CC20=%s.dat"%(cwd,c20))[1:]
#        ic = zip(*ic); ic = np.argmax(ic[0])
#        ax1.plot(c20, 
#                 np.loadtxt("%s/omega.CC20=%s.dat"%(cwd,c20))[ic+1][0],'go')
        
    omegaR = zip(*omegaR)    #omegaIC not working
    ax1.plot(omegaR[0], omegaR[1], 'b-') 
#    ax2.plot(omegaR[0], omegaR[2], 'r-')
      
    for m in model:
        ofile = np.loadtxt("%s/mcst-%s.dat"%(cwd,m))[2]
#        print(mode, m, ofile, ellip[c], ofile + ellip[c])
        ofile = np.append(ofile + ellip[c],
                          np.loadtxt("%s/omega.%s.dat"%(cwd,m))[0])
        ax1.plot(ofile[0], ofile[1], 's', label=m)
#        ax2.plot(ofile[0], ofile[2], 'd')
        
    ax1.axhline(y=SC[c][0], color='r', linestyle='--', label='My Inversion')
#    ax2.axhline(y=SC[c][1], color='r', linestyle='--', label='My Inversion')
    if mode is not "11s00":
        ax1.axhline(y=REM[c][0], color='m', linestyle='--', label='REM')
#        ax2.axhline(y=REM[c][1], color='m', linestyle='--', label='REM')
    ax1.axhline(y=PREM[c][0], color='k', linestyle='--', label='PREM')
#    ax2.axhline(y=PREM[c][1], color='k', linestyle='--', label='PREM')
    
    if mode in ["02s00", "05s00"]:
        ax1.plot(CCcoef[c]+ellip[c], CC[c][0], 'D', label="CC")
#    plt.ylim(PREM[c][0]-0.001)
    ax1.set_title(mode)
    ax1.legend()#loc=1)
    plt.show()
    c += 1