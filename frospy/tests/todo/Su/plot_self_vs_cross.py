# -*- coding: utf-8 -*-
"""
Created on Sat May 19 17:24:35 2018

@author: sujania
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

inv = []
with open('radial_results_selfcross.txt','rb') as source:
    for line in source:
        lines = line.split('\t')
        inv.append(lines)
#inv = inv[1::]

# ----------------------------------------------------
#color = ['b', 'r', 'm', 'g']; 
#width = [-0.15, -0.05, 0.05, 0.15]
color = ['skyblue', 'royalblue', 'navy']; 
width = [-0.05, 0.05, 0.15]

modelname = [r'self-coupling', r'cross-coupling', r'Laske & Masters']
i=0
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(2,1,1)
for model in ['self', 'cross', 'LM']:
    l =[[],[]]
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        l[0].append(float(m[0])); l[1].append(float(m[15]))
    if model=='self' or model=='cross':
        ax.bar(np.array(l[0])+width[i],l[1], 0.1, color=color[i])
    else:
        ax.bar(np.array(l[0])+width[i],l[1], 0.1, color=color[i], yerr=float(m[6]),
               error_kw=dict(ecolor=color[i], lw=2, capsize=4, capthick=2))
    i=i+1
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=3, zorder=0)
plt.text(3.5, -1.5, 'PREM', fontsize=12, backgroundcolor='w', va='bottom', ha='center')
ax.tick_params(axis = 'both', which = 'major', labelsize=12)
ax.yaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylabel('$\delta$f ($\mu$Hz)', fontsize=20)
ax.set_ylim(-4,10)
ax.set_xticklabels([])
legend=[]
for i in range(len(modelname)): 
    legend.append(Line2D([0],[0], color='w', marker='s', markersize=15, 
                       mfc=color[i], label=modelname[i]))#mec=color[i], mfc='w'))
ax.legend(handles=legend, loc='upper left', 
          fontsize=15, frameon=False)

i=0
ax = fig.add_subplot(2,1,2)
for model in ['self', 'cross']:
    l =[[],[]]
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        l[0].append(float(m[0])); l[1].append(float(m[16]))
    ax.bar(np.array(l[0])+width[i],l[1], 0.1, color=color[i])
    i=i+1
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=3, zorder=0)
plt.text(3.5, -150, 'PREM', fontsize=12, backgroundcolor='w', va='bottom', ha='center')
ax.tick_params(axis = 'both', which = 'major', labelsize = 12, zorder=0)
ax.yaxis.set_major_locator(plt.MaxNLocator(7))
ax.set_xlabel('overtone number $n$', fontsize=20)
ax.set_ylabel('$\delta$$Q$',fontsize=20)
ax.set_ylim(-350,1100)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height])
plt.show()
#fig.savefig('self_cross_results', dpi=350) 
