#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 15:02:42 2019

@author: talavera
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

xmin=1000
xmax=4500
ymin=0
ymax=135.75

ifiles = []
ifiles.append("adiabat.in")
ifiles.append("nonadiabat.out")

models = []
for f in ifiles:
    m = np.loadtxt(f)
    m = [ii for ii in zip(*m)]
    models.append(m)

colors = ["r", "b"]    

fig = plt.figure(figsize=(7.5, 12))
ax = fig.add_subplot(111)
plt.gca().invert_yaxis()

for f, m, c in zip(ifiles, models, colors):
    f1 = f.split('.')[0]
    f2 = f.split('.')[1]
    ax.plot(m[1], np.array(m[0])/10000, 
            color=c, label="%sput %s structure" % (f2,f1)) 

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymax,ymin)

ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.tick_params(axis="both", which='both', bottom=True, top=True, right=True,
               labelbottom=True, labeltop=True)
ax.locator_params(axis='x',nbins=13)
ax.locator_params(axis='y',nbins=14)
ax.set_xlabel("temperature (K)")
ax.set_ylabel("pressure (GPa)")
plt.legend()
plt.grid()
plt.savefig('nonadiabat.ps',
              dpi=400,
              pad_inches=0.05,
              )
plt.show()
