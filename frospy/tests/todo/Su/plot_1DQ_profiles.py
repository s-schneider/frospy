#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:26:53 2020

@author: talavera
"""

import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
from os.path import join, basename, splitext
from frospy.util.base import sort_human

dhome = "/net/home/talavera/codes/mineos/DEMO/models/"
                        
ifiles = glob.glob(join(dhome,"prem_ocean.UMM*"))
ifiles = sort_human(ifiles)
ifiles.insert(0,ifiles.pop(1)) # changinng PREM location
ifiles.insert(1,ifiles.pop(3)) # changing QM1 location

Qk = []
Qm = []
name = []
for f in ifiles: 
    n = f.split("/")[-1].split("_")[-2].split("UMM")[1]
    f = open(f, 'r') 
    lines = f.read().splitlines()
    R  = []
    qk = []
    qm = []
    for l in lines[3::]:
        l = l.split()
        R.append(-float(l[0])/1000 + 6371.0)
        qk.append(float(l[4]))
        qm.append(float(l[5]))
    Qk.append(qk)
    Qm.append(qm)
    name.append(n)


colors = ["r", "b", "k", "k", "grey", "grey",  ]
lines = ["-", "-", "-", "--", "-", "--",  ]

#Qkappa & Qkappa
fig = plt.figure()
fig, (ax1, ax, axt) = plt.subplots(1,3,sharey=False,
                             gridspec_kw={'width_ratios': [1,0.7,0.2]},
                             facecolor='w')
fig.set_size_inches(6, 6)

xaxs_split=[10**2,9**6,10**8,10**9,]
ax.spines['right'].set_visible(False)
axt.spines['left'].set_visible(False)
ax.set_xlim(xaxs_split[0], xaxs_split[1])
axt.set_xlim(xaxs_split[2], xaxs_split[3])
axt.tick_params(left=False)

# directly from:
# https://matplotlib.org/examples/pylab_examples/broken_axis.html
d = .01  # how big to make the diagonal lines in ax
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
axt.plot((1-d,1+d),(1-d,1+d),**kwargs) # top/left
ax.plot((1-d,1+d),(-d,+d),**kwargs) # bottom/left

kwargs.update(transform=axt.transAxes, color='k', clip_on=False)
axt.plot((-d,+d),(1-d,1+d),**kwargs) # top/right
ax.plot((-d,+d),(-d,+d),**kwargs) # bottom/right

# plotting        
for n, qk, qm, c, l in zip(name, Qk, Qm, colors, lines):
    ax1.plot(qm, R, label=n, color=c, linestyle=l)
    ax.plot(qk, R, label=n, color=c, linestyle=l)
    axt.plot(qk, R, label=n, color=c, linestyle=l)
    
#Qmu
ax1.legend(frameon=False)
ax1.set_xlabel("$Q_\mu$", fontsize=12)  
ax1.set_ylabel("Depth (km)") 
extraticks=[80, 220, 400, 670, 1470, 2170, 2891, 5149, 5770]
ax1.set_yticks(extraticks) #list(ax1.get_yticks()) + 
ax1.tick_params(axis='both', right=True,which='major', labelsize=10)
ax1.set_ylim(6371, 0)
#ax1.set_xlim(-30, 630)
#ax1.hlines(y=220, xmin=-30, xmax=630, color="lightgray", zorder=0)
#ax1.hlines(y=410, xmin=-30, xmax=630, color="lightgray", zorder=0)
#ax1.hlines(y=660, xmin=-30, xmax=630, color="lightgray", zorder=0)

#Qkappa
ax.set_yticklabels([])
ax.set_xscale('log')    
ax.set_ylim(6371, 0)
ax.set_xlim(5**3, 15**5)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.set_xticks([100, 1000, 60000]) 
ax.set_yticks(extraticks)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel("$Q_\kappa$", fontsize=12)  

axt.set_xscale('log')    
#axt.xaxis.set_major_formatter(NullFormatter())
axt.xaxis.set_minor_formatter(NullFormatter())
axt.yaxis.set_major_formatter(NullFormatter())
axt.yaxis.tick_right()
axt.set_yticks(extraticks) 
#axt.set_title("$\infty$")
axt.set_ylim(6371, 0)

plt.text(11.5**7,   50.0, "...", fontsize=15)
plt.text(11.5**7,  670.0, "...", fontsize=15)
plt.text(11.5**7, 2900.0, "...", fontsize=15)
plt.text(11.5**7, 5149.5, "...", fontsize=15, color="b")
plt.text(10.5**9, 700, "Upper\nMantle", fontsize=12, rotation=-90)
plt.text(10.5**9, 2500, "Lower Mantle", fontsize=12, rotation=-90)
plt.text(10.5**9, 4500, "Outer Core", fontsize=12, rotation=-90)
plt.text(10.5**9, 6300, "Inner Core", fontsize=12, rotation=-90)
#plt.suptitle("1D attenuation models",y=0.92)
plt.suptitle("Our measurements",)#fontweight="bold")
plt.show() 
fig.savefig('1DQ.png',orientation='landscape', 
            dpi=400, bbox_inches='tight', pad_inches=0.02)
