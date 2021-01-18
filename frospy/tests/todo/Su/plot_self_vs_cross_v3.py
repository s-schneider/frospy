# -*- coding: utf-8 -*-
"""
Created on Sat May 19 17:24:35 2018

@author: sujania
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

inv = []
home = '/Users/sujania/Documents/'
results = '%s/PhD/NotesMAC/digitalizations/Lynthgoe15/radialshift.txt' % home
with open(results,'rb') as source:
    for line in source:
        lines = line.split('\t')
        inv.append(lines)
#inv = inv[1::]

# ----------------------------------------------------
color = ['skyblue', 'hotpink','b', 'r']; 
width = [-0.15, -0.05, 0.05, 0.15]
modelname = [r'Laske & Masters, 2001: self-coupling', 
             r'Laske & Masters, 2001: cross-coupling',
             r'My Inversion: self-coupling', 
             r'My Inversion: cross-coupling']
i=0
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(2,1,1)
for model in ['LM-SC', 'LM-CC', 'SC', 'CC']:
    l =[[],[]]
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        l[0].append(float(m[0])); l[1].append(float(m[7]))
        if model=='SC' or model=='CC':
            ax.errorbar(int(m[0])+width[i], float(m[7]), fmt='*',
                        ms=15, color=color[i], capsize=5, elinewidth=2.5)
                        #, fillstyle='none'
        else:
            ax.errorbar(int(m[0])+width[i], float(m[7]), yerr=float(m[6]), 
                        fmt='s', ms=9, color=color[i], capsize=5, 
                        elinewidth=2.5) #, fillstyle='none')
    i=i+1
print np.linspace(1, len(mode), len(mode))
for j in np.linspace(1, len(mode), len(mode)):
    plt.axvline(x=j+0.5, color='lightgray', 
                linestyle='--',  linewidth=1, zorder=0)
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=5, zorder=0)
plt.text(10, 0, 'PREM', fontsize=19, backgroundcolor='w', 
         va='center', ha='center', zorder=1)
ax.tick_params(axis = 'both', which = 'major', labelsize=18, zorder=0)
modename = ['${}_{%s}S_0$'%int(ll) for ll in l[0]]
ax.set_xticks(np.linspace(1, len(mode), len(mode)))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=25)
xlim = ax.get_xlim()
ax.yaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylabel('$\delta$f ($\mu$Hz)', fontsize=25)
ax.set_ylim(-3,8.5)

i=0
color = ['b', 'r']
width = [-0.15, -0.05]
#modelname = [r'self-coupling', r'cross-coupling']
ax = fig.add_subplot(2,1,2)
for model in ['SC', 'CC']:
    l =[[],[]]
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        l[0].append(float(m[0])); l[1].append(float(m[11]))
        ax.errorbar(int(m[0])+width[i], float(m[11]), fmt='*',
                    ms=15, color=color[i], capsize=5, elinewidth=2.5)
        #, fillstyle='none')
    i=i+1
print np.linspace(1, len(mode), len(mode))
for j in np.linspace(1, len(mode), len(mode)):
    plt.axvline(x=j+0.5, color='lightgray', linestyle='--',  
                linewidth=1, zorder=0)
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=3, zorder=0)
plt.text(10, 0, 'PREM', fontsize=19, backgroundcolor='w', 
         va='center', ha='center')
ax.tick_params(axis = 'both', which = 'major', labelsize = 18, zorder=0)
ax.yaxis.set_major_locator(plt.MaxNLocator(7))
ax.set_ylabel('$\delta$$Q$',fontsize=20)
ax.set_ylim(-350,500)
modename = ['${}_{%s}S_0$'%int(ll) for ll in l[0]]
ax.set_xticks(np.linspace(1, len(mode), len(mode)))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=25)
xlim = ax.get_xlim()
ax.yaxis.set_major_locator(plt.MaxNLocator(4))

legend=[]
color = ['b', 'r', 'skyblue', 'hotpink']; 
for i in range(len(modelname)): 
    if not(modelname[i]=='My Inversion: self-coupling' or 
           modelname[i]=='My Inversion: cross-coupling'): 
        legend.append(Line2D([0],[0], color='w', marker='s', markersize=15, 
                      mfc=color[i], label=modelname[i]))
                      #mec=color[i], mfc='w'))
    else:
        legend.append(Line2D([0],[0], color='w', marker='*', markersize=25, 
                      mfc=color[i], label=modelname[i]))
                      #mec=color[i], mfc='w'))
lgd = ax.legend(handles=legend, loc='lower center', 
                fontsize=14, frameon=False, ncol=2,  
                bbox_to_anchor=(0.45, -0.5))

plt.show()
fig.savefig('self_cross_results',  bbox_extra_artists=(lgd,), 
            bbox_inches='tight', dpi=350) 