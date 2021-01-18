#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:43:21 2019

@author: talavera
"""
import numpy as np
import matplotlib.pyplot as plt
from frospy.core.modes import read as read_modes
from frospy.core.splittingfunc.plot import sens_kernel

allmodes = read_modes()

mtype = 'l2'
if mtype == 'RM':
    ## Radial modes
    #mrange = range(1,12)
    #mrange = [0, 1, 2, 3, 4, 5, 6, 8, 9]
    mrange = [0, 1, 2, 0, 3, 0, 0, 4, 0, 5, 6, 7, 8, 9, 11]
    modenames = ['${}_{%s}S_0$'%int(ll) for ll in mrange]
    modes = ['%02dS00'%int(ll) for ll in mrange]
#    width = np.ones(len(mrange))*0.5
    width = []
    width.extend(np.ones(3)*0.5)
    width.extend(np.ones(1)*0.35)
    width.extend(np.ones(1)*0.5)
    width.extend(np.ones(1)*0.25)
    width.extend(np.ones(1)*0.35)
    width.extend(np.ones(1)*0.5)
    width.extend(np.ones(1)*0.35)
    width.extend(np.ones(6)*0.5)
    gridspec_kw = {'width_ratios': width, 'wspace':0.}
    fig, axes = plt.subplots(nrows=1, ncols=len(mrange), 
                             sharey=True, sharex=False,
                             gridspec_kw=gridspec_kw, )

elif mtype == 'l2':
    ##IC modes
    #mrange = [4, 7, 8, 9, 10, 11, 13, 27]
    #mrange = [2, 4, 7, 8, 9, 10, 11, 13, 16, 18, 20, 22, 25, 27]
    mrange = [2, 4, 7, 0, 8, 9, 0, 10, 11, 0, 13, 16, 18, 20, 22, 27]
    modenames = ['${}_{%s}S_2$'%int(ll) for ll in mrange]
    modes = ['%02dS02'%int(ll) for ll in mrange]

#    width = np.ones(len(mrange))*0.5
    width = []
    width.extend(np.ones(3)*0.5)
    width.extend(np.ones(1)*0.1)
    width.extend(np.ones(2)*0.5)
    width.extend(np.ones(1)*0.1)
    width.extend(np.ones(2)*0.5)
    width.extend(np.ones(1)*0.1)
    width.extend(np.ones(6)*0.5)
    gridspec_kw = {'width_ratios': width, 'wspace':0.}
    fig, axes = plt.subplots(nrows=1, ncols=len(modes), 
                             sharey=True, sharex=False,
                             gridspec_kw=gridspec_kw, )
fig.set_size_inches(3.74, 1.4)

legend_show = False       
ticks = True # -intmax*0.4
text = True
fontsize=4
i = 0
for ax, mode, name in zip(axes.flat, modes, modenames):   
    m = allmodes.select(mode)
    if text:
        if m[0].n == 0 and m[0].l == 0: 
    #        ax.text(-18, 4000, "Depth (km)", fontsize=10, rotation=90)
    #        legend_show = True
    #        ticks = False
            ax.text(-7., 5700, "TZ",  fontsize=fontsize, va="bottom")
            ax.text(-7., 3480, "CMB", fontsize=fontsize, va="bottom")
            ax.text(-7., 1220, "ICB", fontsize=fontsize, va="bottom")
            ax.text(-16, 1200, "Radial Modes", fontsize=8, rotation=90)
            ax.text(-16, 6600, "a)", fontsize=8)
#            ax.spines['left'].set_linewidth(1.5)
            text=False
        elif m[0].n == 2 and m[0].l == 2: 
    #        ax.text(-20, 4000, "Depth (km)", fontsize=10, rotation=90)
    #        legend_show = True
            ax.text(-17, 5700, "TZ",  fontsize=fontsize, va="bottom")
            ax.text(-17, 3480, "CMB", fontsize=fontsize, va="bottom")
            ax.text(-17, 1220, "ICB", fontsize=fontsize, va="bottom")
            ax.text(-36, 1000, "$l=2$ S. Modes", fontsize=8, rotation=90)
#            ax.text(-40, 1000, "Splitting function coefficients of IC sensitive modes", fontsize=10, weight="bold",rotation=90)
            ax.text(-36, 6600, "b)", fontsize=8)
#            ax.spines['left'].set_linewidth(1.5)
            text=False
            
#    for axis in ['top','bottom']:
#        ax.spines[axis].set_linewidth(1.5)
#    print(mode,m[0].l,i)
    ii = [0,1,2,4,7,9,10,11,12,13,14]
    jj = [0,1,2,4,5,7,8,10,11,12,13,14,15]
    if (m[0].l == 2 and i in jj) or (m[0].l == 0 and i in ii):
#        print(m[0].l,i)
        sens_kernel(m[0], title=name, ax=ax, legend_show=legend_show,
                    ticks = ticks, fontsize=fontsize, linewidth=0.6)
    elif m[0].l == 0 and i in [6]:
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
    elif m[0].l == 0 and i in [5]:
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
    if m[0].l == 2 and i in [5,8]:
        print(mode)
        ax.spines['left'].set_visible(False)

    ax.set_xticklabels([])
    ax.set_xticks([])
    legend_show = False       
    ticks = True
    i = i + 1
#ax.spines['right'].set_linewidth(1.5)
ax.axes.get_xaxis().set_ticks([])
ax.axes.get_yaxis().set_ticks([])

#fig.savefig('%s_senskernels.png'%mtype,orientation='landscape', 
#            dpi=400, bbox_inches='tight', pad_inches=0.02)

## one mode
#from frospy.core.modes import read as read_modes
#from frospy.core.splittingfunc.plot import sens_kernel
#
#name = '16S7'
#m = read_modes(modenames=name)
#sens_kernel(m[0], title=name)