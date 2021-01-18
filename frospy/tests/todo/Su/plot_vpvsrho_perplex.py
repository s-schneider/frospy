#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 10:40:45 2019

@author: talavera
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

# Plots Vp, Vs, rho as a function of depth

f=open('perplex_setup.dat')
lines=f.read().splitlines()

xmin=float(lines[20].split()[0])
xmax=float(lines[20].split()[1])
ymin=float(lines[22].split()[0])
ymax=float(lines[22].split()[1])

ifiles = []
ifiles.append("prem.dat")
ifiles.append("ak135")
ifiles.append(lines[25].split()[0])
ifiles.append(lines[26].split()[0])
ifiles.append(lines[1]+".vpvsrho")

models = []
for f in ifiles:
    m = np.loadtxt(f)
    m = [ii for ii in zip(*m)]
    models.append(m)

colors = ["r", "r", "g", "k", "b"]
fig = plt.figure(figsize=(7.5, 12))
ax = fig.add_subplot(111)
plt.gca().invert_yaxis()

for f, m, c in zip(ifiles, models, colors):
    f = f.split('.')[0]
    if f == "prem":
        ax.plot(m[2], m[0], color=c, label=f) # Vp
        ax.plot(m[3], m[0], color=c) # Vs
        ax.plot(m[1], m[0], color=c) # density
    elif f == "ak135":
        ax.plot(m[2], m[0], color=c, linestyle='--', label=f) # Vp
        ax.plot(m[3], m[0], color=c, linestyle='--') # Vs
        ax.plot(m[1], m[0], color=c, linestyle='--') # density
    else:
        ax.plot(m[2], m[5], color=c, label=f) # Vp
        ax.plot(m[3], m[5], color=c) # Vs
        ax.plot(np.array(m[4])/1000, m[5], color=c) # density

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymax,ymin)
ax.text(x=4, y=2400, s="density \n (g/cc)", color="r")
ax.text(x=7, y=1400, s="$v_s$ \n (km/s)", color="r")
ax.text(x=12.5, y=1300, s="$v_p$ \n (km/s)", color="r")
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))
ax.tick_params(axis="both", which='both', bottom=True, top=True, right=True,
               labelbottom=True, labeltop=True)
ax.locator_params(axis='x',nbins=xmax-xmin)
ax.locator_params(axis='y',nbins=20)
ax.set_ylabel("depth (km)")
plt.legend()           
plt.grid()
plt.savefig('%s.ps'%ifiles[-1].split('.')[0],
              dpi=400,
              pad_inches=0.05,
              )
plt.show()
