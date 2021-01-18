# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 16:30:40 2018

@author: sujania
"""

# -*- coding: utf-8 -*-
"""
Created on Sat May 19 17:24:35 2018

@author: sujania
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from colour import Color

inv = []
home = '/Users/sujania/Documents/'
#results = '%s/PhD/NotesMAC/digitalizations/Lynthgoe15/radialshift.txt' % home
results = '%s/PhD/NotesUU/RadialModes/results_AD/radialshift.txt' % home
with open(results,'rb') as source:
    for line in source:
        lines = line.split('\t')
        inv.append(lines)
#inv = inv[1::]

# ----------------------------------------------------
reds = plt.cm.get_cmap("Greens", 10)
blues = plt.cm.get_cmap("Blues", 10)


color = [blues(7), blues(5), blues(2), 
         reds(2), reds(4), reds(6), 
#         'r', 'dodgerblue', 'b', 'm', 'b']
         'r', 'b', 'm', 'b']

modelname = [r'$\phi=1.20$', 
             r'$\phi=1.10$', 
             r'$\phi=1.04$', 
             r'$\phi=0.96$', 
             r'$\phi=0.90$', 
             r'$\phi=0.80$',
             r'My Inversion: SC', #: self-coupling', 
             r'My Inversion: CC', 
#             r'My Inversion: CC with ellip', 
#             r'My Inversion: CC with structure', 
#             r'Widmer et al, 1991', 
#             r'Masters & Widmer (REM)', 
#             r'Durek & Ekstrom, 1995',
#             r'He & Tromp, 1996', 
             ]
i=0
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(1,1,1)
for model in ['phi=1.20', 'phi=1.10', 'phi=1.04', 
              'phi=0.96', 'phi=0.90', 'phi=0.80',
              'SC', 'CC']:#, 'CC']:#,'QM1', 'REM', ]:
    l =[[],[]]
    mode = [m for m in inv if m[3]==model]
    for m in mode:
        l[0].append(float(m[0])); l[1].append(float(m[7]))
        if model in ['SC','SCC','CC']:
            ax.errorbar(int(m[0]), float(m[7]), fmt='*',
                        ms=20, color=color[i], capsize=5, elinewidth=2.5)

        elif model=='REM' or model=='HT' or model=='QM1' or model=='DE':
            ax.errorbar(int(m[0]), float(m[7]), yerr=float(m[6]), 
                        fmt='s', ms=11, color=color[i], capsize=5, 
                        elinewidth=2.5) 
        else:
            ax.errorbar(int(m[0]), float(m[7]), fmt='o', ms=7, 
                        color=color[i], capsize=5, 
                        elinewidth=2.5)
    i=i+1
for j in np.linspace(1, len(mode), len(mode)):
    plt.axvline(x=j+0.5, color='lightgray', linestyle='--', 
                linewidth=1, zorder=0)
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=3, zorder=0)
plt.text(7.5, 0, 'PREM', fontsize=16, backgroundcolor='w', 
         va='center', ha='center')
ax.tick_params(axis = 'both', which = 'major', labelsize=18)
ax.yaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylabel('$\delta$f ($\mu$Hz)', fontsize=20)
ax.set_ylim(-20,15)
#modename = ['${}_{%s}S_0$'%str(j) for j in range(1, len(mode)+1)]
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(1,10)]
ax.set_xticks(np.linspace(1, len(modename)+1, len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=25)
ax.set_xlim(0.5,9.5)

legend=[]
for i in range(len(modelname)): 
    if modelname[i] in ['He & Tromp, 1996', 'Masters & Widmer (REM)']:
        legend.append(Line2D([0],[0], color='w', marker='s', markersize=15, 
                      mfc=color[i], label=modelname[i]))
    elif modelname[i] in ['My Inversion: SC', 'My Inversion: CC']: # with ellip', 
#                       'My Inversion: CC with structure']: 
        legend.append(Line2D([0],[0], color='w', marker='*', markersize=25, 
                      mfc=color[i], label=modelname[i]))
    else:
        legend.append(Line2D([0],[0], color='w', marker='o', markersize=15, 
                      mfc=color[i], label=modelname[i]))
lgd = ax.legend(handles=legend, loc='lower center', 
                fontsize=14, frameon=False, ncol=3,  
                bbox_to_anchor=(0.45, -0.35))

plt.show()
fig.savefig('fshift_radial_ani',  bbox_extra_artists=(lgd,), 
            bbox_inches='tight', dpi=350) 