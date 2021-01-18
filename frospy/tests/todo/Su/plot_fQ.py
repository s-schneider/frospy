#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:27:41 2019

@author: talavera
"""

#--------------------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
from frospy.core.modes import read as read_modes
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

mpl.rcParams['axes.linewidth'] = 1.8 #set the value globally
#plt.rc('text', usetex=True)

modes_all = read_modes()
mode = "04s00"
cwd = '/net/home/talavera/eejit/splitting/%s/synthetics/IC/QCC'%mode

fig = plt.figure(figsize=(3.74,4.52))

ax = fig.add_subplot(121)
m1 = 'cstC'
m2 = 'PREMVs'
m3 = 'highVsSelf'
if mode == "04s00":
    reverse=True
else:
    reverse=False

omega2 = np.loadtxt('%s/omega.%s.dat'%(cwd,m2))
omega2 =sorted(omega2, key=lambda x: x[1], reverse=reverse); 
omega2 = zip(*omega2)
p1 = ax.scatter(omega2[0][0:5], omega2[1][0:5], linewidths=0.8,
                color='k', edgecolors='k', marker='D', s=12, zorder=50,
                label='PREM') #'dimgray',
p2 = ax.scatter(omega2[0][5::], omega2[1][5::], linewidths=0.8,
                color='w', edgecolors='k', marker='D', s=12, zorder=50)

omega3 = np.loadtxt('%s/omega.%s.dat'%(cwd,m3))
omega3 =sorted(omega3, key=lambda x: x[1], reverse=reverse); 
omega3 = zip(*omega3)
p3 = ax.scatter(omega3[0][0:5], omega3[1][0:5], linewidths=0.8,
                color='grey', edgecolors='grey', marker='o', s=20,
                label='$+2 \%$ $v_s$')
p4 = ax.scatter(omega3[0][5::], omega3[1][5::], linewidths=0.8,
                color='w', edgecolors='grey', marker='o', s=20)

omega1 = np.loadtxt('%s/omega.%s.dat'%(cwd,m1))
omega1 =sorted(omega1, key=lambda x: x[1], reverse=reverse); 
omega1 = zip(*omega1)
p5 = ax.scatter(omega1[0][0:5], omega1[1][0:5], linewidths=0.8,
                color='r', edgecolors='r', marker='o', s=20, zorder=100,
                label='Our Inversion')
p6 = ax.scatter(omega1[0][5::], omega1[1][5::], linewidths=0.8,
                color='w', edgecolors='r', marker='o', s=20, zorder=100)

modesin = open('%s/modes.in'%cwd).read().splitlines()
for i, mode in enumerate(modesin):
    modesin[i] = mode.split()
    
title = []
for mode in modesin[1:int(modesin[0][0])+1]:
    mode = modes_all.select(name=''.join(mode))[0]
    label = '${}_{%s}%s_{%s}$' % (mode.n, mode.type, mode.l)
    title.append(label)
#    ax.text(mode.freq, mode.Q*1.15, label, fontsize=20)
    if mode.n == 10:
        y = mode.Q*1.5
        x = mode.freq*0.999
        ax.set_ylim([-20,870])
        num = 'a)'

        l = ax.legend([(p1, p2), (p3, p4), (p5, p6)], 
                      ["PREM", "$+2\%$ $v_s$", "Our inversion"],
                      handler_map={tuple: HandlerTuple(ndivide=None)},
                      loc='upper left', bbox_to_anchor=(-0.01, 0.98), 
                      columnspacing=0,handletextpad=0.25,
                      labelspacing=0.1,handlelength=3,
                      ncol=1, fontsize=8, frameon=False)        
        
        xl = ax.get_xlim()[0]
        ax.text(xl*1.0025, 795, title[0], fontsize=7, zorder=600)
        
#        ax.legend(loc='upper left', bbox_to_anchor=(-0.05, 1.03), 
#                  columnspacing=0,handletextpad=0.,
#                  labelspacing=0.1,
#                  ncol=1, fontsize=15, frameon=False) 
#        xl = ax.get_xlim()[0]*1.003
#        ax.scatter(xl,780, linewidths=2,
#                   color='w', edgecolors='grey', marker='D', s=30, )
#        ax.scatter(xl,700, linewidths=2,
#                   color='w', edgecolors='k', marker='o', s=50, )
#        ax.scatter(xl,600, linewidths=2,
#                   color='w', edgecolors='r', marker='o', s=50, )
    if mode.n == 11:
        y = mode.Q*0.11
        x = mode.freq*0.998
        ax.set_title(#'No 1D $Q$ coupling',
                     'No 1D Q coupling \n $\mathrm{%s}$ %s-%s'%(num,title[0],title[1]), 
                     fontsize=10, weight="bold" )
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
        ax.text(xl*1.0085, 795, title[1], fontsize=7, zorder=600)  
    if mode.n == 8:
        y = mode.Q*0.55
        x = mode.freq*0.997
        ax.set_ylim([150,600])
        num = 'c)'
        l = ax.legend([(p1, p2), (p3, p4), (p5, p6)], 
                      ["", "", ""],
                      handler_map={tuple: HandlerTuple(ndivide=None)},
                      loc='upper left', bbox_to_anchor=(-0.01, 0.98), 
                      columnspacing=0,handletextpad=0.25,
                      labelspacing=0.1,handlelength=3,
                      ncol=1, fontsize=8, frameon=False)        
        
    if mode.n == 9:
        y = mode.Q*1.2
        x = mode.freq*0.998
        ax.set_title('%s %s-%s'%(num,title[0],title[1]), 
                     fontsize=10)
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.005)) 
    ax.text(x, y, label, fontsize=10, zorder=600)
ax.set_xlabel("f (mHz)", fontsize=8)
ax.set_ylabel("Q", fontsize=8)    
ax.tick_params(axis='both', which='major', labelsize=8)
ax.xaxis.set_major_locator(plt.MaxNLocator(3))
#ax.set_title('a) No 1D coupling', fontsize=18)


ax2 = fig.add_subplot(122)
m1 = 'cstCatt'
m2 = 'PREMVsCatt'
m3 = 'highVsCatt'

omega2 = np.loadtxt('%s/omega.%s.dat'%(cwd,m2))
omega2 =sorted(omega2, key=lambda x: x[1], reverse=reverse); 
omega2 = zip(*omega2)
ax2.scatter(omega2[0][0:5], omega2[1][0:5], linewidths=0.8,
            color='k', edgecolors='k', marker='D', s=12, zorder=50)
ax2.scatter(omega2[0][5::], omega2[1][5::], linewidths=0.8,
            color='w', edgecolors='k', marker='D', s=12, zorder=50) 
            #'dimgray',

omega3 = np.loadtxt('%s/omega.%s.dat'%(cwd,m3))
omega3 =sorted(omega3, key=lambda x: x[1], reverse=reverse); 
omega3 = zip(*omega3)
ax2.scatter(omega3[0][0:5], omega3[1][0:5], linewidths=0.8,
            color='grey', edgecolors='grey', marker='o', s=20)
ax2.scatter(omega3[0][5::], omega3[1][5::], linewidths=0.8,
            color='w', edgecolors='grey', marker='o', s=20)

omega1 = np.loadtxt('%s/omega.%s.dat'%(cwd,m1))
omega1 =sorted(omega1, key=lambda x: x[1], reverse=reverse); 
omega1 = zip(*omega1)
ax2.scatter(omega1[0][0:5], omega1[1][0:5], linewidths=0.8,
            color='r', edgecolors='r', marker='o', s=20, zorder=100)
ax2.scatter(omega1[0][5::], omega1[1][5::], linewidths=0.8,
            color='w', edgecolors='r', marker='o', s=20, zorder=100)

# adjusting plot so both have the ame limits and aspect
xlim = ax2.get_xlim()
ylim = ax.get_ylim()

ax.set_xlim(xlim)
aspect = (ax.get_xlim()[1] - ax.get_xlim()[0]) /   \
         (ax.get_ylim()[1] - ax.get_ylim()[0])
ax.set_aspect(adjustable='box', aspect=aspect)

ax2.set_ylim(ylim)
for mode in modesin[1:int(modesin[0][0])+1]:
    mode = modes_all.select(name=''.join(mode))[0]
    label = '${}_{%s}%s_{%s}$' % (mode.n, mode.type, mode.l)
#    ax2.text(mode.freq, mode.Q*1.15, label, fontsize=20)
    if mode.n == 10:
        y = mode.Q*3.2
        x = mode.freq*1.005
        num = 'b)'
    if mode.n == 11:
        y = 5#mode.Q*0
        x = mode.freq*0.999
        ax2.set_title(#'1D $Q$ coupling',
                      '1D Q coupling \n $\mathrm{%s}$ %s-%s'%(num,title[0],title[1]), 
                      fontsize=10, weight="bold")#, loc="left")
        ax2.yaxis.set_minor_locator(ticker.MultipleLocator(50))
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
    if mode.n == 8:
        y = mode.Q*0.55
        x = mode.freq*0.995
        num = 'd)'
        ax2.set_title('%s %s-%s'%(num,title[0],title[1]), 
                      fontsize=10) 
        ax2.yaxis.set_minor_locator(ticker.MultipleLocator(25))
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
        xl = ax.get_xlim()[0]
        ax.text(xl*1.001, 565, title[0], fontsize=7, zorder=600)
    if mode.n == 9:
        y = mode.Q*1.35
        x = mode.freq
        ax.text(xl*1.004, 565, title[1], fontsize=7, zorder=600) 
    ax2.text(x, y, label, fontsize=10, zorder=200)
aspect = (ax2.get_xlim()[1] - ax2.get_xlim()[0]) /   \
         (ax2.get_ylim()[1] - ax2.get_ylim()[0])
ax2.set_aspect(adjustable='box', aspect=aspect)
ax2.set_xlabel("f (mHz)", fontsize=8)
ax2.tick_params(axis='both', which='major', labelsize=8)
ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
ax2.tick_params(axis='y',labelleft='off') 
ax2.yaxis.set_ticks_position('none') 

plt.subplots_adjust(wspace=0.05)
#plt.savefig('fQ_%s.png'%mode.name,
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.05,
#              )
plt.show()