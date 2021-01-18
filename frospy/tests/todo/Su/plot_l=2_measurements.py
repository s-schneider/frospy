#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 14:06:10 2019

@author: talavera
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May 20 10:15:34 2018

@author: sujania
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

inv = []
home = '/Users/sujania/Documents/'
results = '%s/PhD/NotesUU/RadialModes/results_AD/l2shift.txt' % home
with open(results,'rb') as source:
    for line in source:
        lines = line.split('\t')
        inv.append(lines)

# ----------------------------------------------------
mnum = 6
l2modes = [4,7,9,10,13,27]
color = ['gray', 'gray', 'gray', 'gray', 'r','b'];
width = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35]
marker = ['v', 'o', 's', 'D'];
modelname = [r'Masters & Widmer, 1995',
             r'He & Tromp, 1996',  
             r'Masters & Widmer (REM)',
             r'Deuss et al, 2013',
             r'My Inversion CC: ellip',
             r'My Inversion CC',
             ]

i=0
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(2,1,1)
for model in ['MW95', 'HT', 'REM', 'AD', 'SCC', 'CC']:
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        index = l2modes.index(int(m[0]))+1
        if model not in ['SCC', 'CC']: 
            ax.errorbar(index+width[i], float(m[7]), yerr=float(m[6]), 
                        marker=marker[i], markersize=9, color=color[i], 
                        capsize=5, elinewidth=2.5)
        else:
            ax.errorbar(index+width[i], float(m[7]), marker='*',
                        markersize=15, color=color[i], capsize=5, 
                        elinewidth=2.5)
    i=i+1
    
modename = ['${}_{%s}S_2$'%int(ll) for ll in l2modes]

for j in np.linspace(1, len(modename), len(modename)):
    plt.axvline(x=j+0.5, color='lightgray', 
                linestyle='--',  linewidth=1, zorder=0)
plt.axhline(y=0, color='lightgray', linestyle='-',  linewidth=5, zorder=0)
plt.text(3.5, 0, 'PREM', fontsize=19, backgroundcolor='w', 
         va='center', ha='center', zorder=1)
ax.tick_params(axis = 'both', which = 'major', labelsize=18, zorder=0)
ax.set_xticks(np.linspace(1, len(modename), len(modename)))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=25)
xlim = ax.get_xlim()
ax.yaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylabel('$\delta$f ($\mu$Hz)', fontsize=25)
#ax.set_ylim(-2,6)
ax.set_xlim(0.5,mnum+0.5)

i = 0
ax = fig.add_subplot(2,1,2)
ax.set_xticks(np.linspace(1, len(modename), len(modename)))

for model in ['MW95', 'HT', 'REM', 'AD', 'SCC', 'CC']:
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        index = l2modes.index(int(m[0]))+1
        if model not in ['SCC', 'CC']: 
            ax.errorbar(index+width[i], float(m[11]), yerr=float(m[10]), 
                        marker=marker[i], markersize=9, color=color[i], 
                        capsize=5, elinewidth=2.5)
        else:
            ax.errorbar(index+width[i], float(m[11]), marker='*',
                        markersize=15, color=color[i], 
                        capsize=5, elinewidth=2.5)
    i=i+1
for j in np.linspace(1, len(modename), len(modename)):
    plt.axvline(x=j+0.5, color='lightgray', 
                linestyle='--',  linewidth=1, zorder=0)
plt.axhline(y=0, color='lightgray', linestyle='-',  linewidth=5, zorder=0)
plt.text(3.5, 0, 'PREM', fontsize=19, backgroundcolor='w', 
         va='center', ha='center', zorder=1)
ax.tick_params(axis = 'both', which = 'major', labelsize = 18, zorder=0)
ax.set_xticklabels(modename, rotation='horizontal', fontsize=25)
ax.set_xlim(xlim)
ax.yaxis.set_major_locator(plt.MaxNLocator(7))
ax.set_ylabel('$\delta$$Q$',fontsize=25)
#ax.set_ylim(-100,180)
ax.set_xlim(0.5,mnum+0.5)

legend=[]; i=0
for model in modelname: 
    if model not in ['My Inversion CC: ellip','My Inversion CC']: 
        legend.append(Line2D([0],[0], color='w', marker=marker[i], 
                      markersize=15, mfc=color[i], label=modelname[i]))
    else:
        legend.append(Line2D([0],[0], color='w', marker='*', markersize=20, 
                      mfc=color[i], label=modelname[i]))
    i = i + 1
lgd = ax.legend(handles=legend, loc='lower center', 
                fontsize=14, frameon=False, ncol=2,  
                bbox_to_anchor=(0.5, -0.75))

plt.show()
plt.tight_layout()
#fig.savefig('df_dQ_results',  bbox_extra_artists=(lgd,), 
#            bbox_inches='tight', dpi=350) 
