#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 18:08:16 2019

@author: talavera
"""
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib as mpl
import matplotlib
import numpy as np

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
mpl.rcParams['axes.linewidth'] = 1.2 #set the value globally
blues = plt.cm.get_cmap("Blues", 10)
color = blues(3)
xplot = "R" # vs, vp, rho, R
yplot ="c00" # c00, freq, Q, 

fig = plt.figure(figsize=(1.5,1.5))
ax = fig.add_subplot(111)
ax.tick_params(axis='both', which='major', labelsize=6)

if xplot == "vs":
    # 10s2 changing vs
    m1 =[ 
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [111.9, 55.9, 28.0, 14.0, -14.0, -28.0, -55.9, -111.9], # c00
        [720, 573, 424, 327, 131, 103, 90, 84],  # Q
        [4040.7, 4039.58, 4038.52, 4037.79, 4032.3, 4023.15, 4003.6, 3962.84], # f
        ]
    param_m1 = 0.7 # match with measurment
    
    # 10s2-11s2 changing vs
    cc =[ 
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [-142.7, -71.3, -35.7, -17.8, 17.8, 35.7, 71.3, 142.7], # c00
        ]
    param_cc = 0.55 # match with measurment

    # 11s2 changing vs
    m2 =[  
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [182.9, 91.5, 45.7, 22.9, -22.9, -45.7, -91.5, -182.9], # c00
        [91, 92, 96, 103, 199, 308, 496, 697], # Q
        [4132.58, 4092.74, 4073.17, 4063.55, 4048.24, 4046.96, 4045.55, 4044.07], # f
        ]
    param_m2 = 0.45 # match with measurment
    
elif xplot == "vp":
    # 10s2 changing vp
    m1 =[ 
        [10, 5, 2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2, -5, -10], # %increase
        [105.5, 52.7, 21.1, 10.5, 5.3, 2.6, -2.6, -5.3, -10.5, -21.1, -52.7, -105.5], # c00
        [86, 94, 218, 221, 220, 220, 221, 221, 222, 224, 230, 237], # Q
        [4055.57, 4051.87, 4045.15, 4040.8, 4038.75, 4037.78, 4035.91, 4034.99, 4033.19, 4029.71, 4019.74, 4003.72], # f
        ]
    param_m1 = 3.8 # match with measurment
    
    # 10s2-11s2 changing vp
    cc =[ 
        [10, 5, 2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2, -5, -10], # %increase
        [17.9, 8.9, 3.6, 1.8, 0.9, 0.4, -0.4, -0.9, -1.8, -3.6, -8.9, -17.9], # c00
        ]
    param_cc = 3.5 # match with measurment

    # 11s2 changing vp
    m2 =[ 
        [10, 5, 2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2, -5, -10], # %increase
        [3.1, 1.6, 0.6, 0.3, 0.2, 0.1, -0.1, -0.2, -0.3, -0.6, -1.6, -3.1], # c00
        [849, 461, 136, 123, 122, 122, 121, 121, 121, 120, 118, 115], # Q
        [4065.90, 4054.36, 4051.93, 4053.21, 4053.73, 4053.94, 4054.28, 4054.43, 4054.69, 4055.11, 4055.87, 4056.48], # f
        ]
    param_m2 = 3.5 # match with measurment

elif xplot == "rho":
    # 10s2 changing rho
    m1 =[
        [40, 20, 10, 5, -5, -10, -20, -40], # %increase
        [12.8, 6.4, 3.2, 1.6, -1.6, -3.2, -6.4, -12.8], # c00
        [669, 497, 366, 294, 155, 121, 100, 90], # Q
        [4022.26, 4030.19, 4033.76, 4035.41, 4036.72, 4034.00, 4026.46, 4010.48], # f
        ]
    param_m1 = 32 # match with measurment

    # 10s2-11s2 changing rho
    cc =[
        [40, 20, 10, 5, -5, -10, -20, -40], # %increase
        [-92.2, -46.1, -23.0, -11.5, 11.5, 23.0, 46.1, 92.2], # c00
        ]
    param_cc = 20 # match with measurment
    
    # 11s2 changing rho
    m2 =[
        [40, 20, 10, 5, -5, -10, -20, -40], # %increase
        [43.9, 22.0, 11.0, 5.5, -5.5, -11.0, -22.0, -43.9], # c00
        [89, 93, 99, 107, 162, 222, 342, 526], # Q
        [4084.50, 4068.72, 4061.18, 4057.54, 4052.23, 4052.94, 4056.45, 4064.29], # f
        ]
    param_m2 = 35 # match with measurment
    
elif xplot == "R":
    # 10s2 changing R
    m1 =[ 
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [-122.1, -61.0, -30.5, -15.3, 15.3, 30.5, 61.0, 122.1], # c00
        [85, 91, 104, 131, 338, 447, 614, 763],  # Q
        [3960.81, 4002.32, 4022.29, 4031.71, 4038.01, 4038.82, 4039.91, 4040.85], # f
        ]
    param_m1 = -0.7 # match with measurment

    # 10s2-11s2 changing R
    cc =[ 
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [148.0, 74.0, 37.0, 18.5, -18.5, -37.0, -74.0, -148.0], # c00
        ]
    param_cc = -0.55 # match with measurment
    
    # 11s2 changing R
    m2 =[
        [2, 1, 0.5, 0.25, -0.25, -0.5, -1, -2], # %increase
        [-174.8, -87.4, -43.7, -21.9, 21.9, 43.7, 87.4, 174.8], # c00
        [652, 467, 298, 198, 102, 95, 91, 90], # Q
        [4045.50, 4046.54, 4047.67, 4048.76, 4063.40, 4073.01, 4092.69, 4132.99], # f
        ]
    param_m2 = -0.45 # match with measurment
       
if yplot == "c00":
    y1  = m1[1]
    y2  = m2[1]
    ycc = cc[1]
    
    # my 10s2 inversion
    obs_m1 = 40.01  
    obs_m1_err = np.array([10, 10]) # obs_m1_err = np.array([-0.694, 3.376])
    obs_m1_err = obs_m1_err.reshape((2, 1))
    
    # my 10s2-11s2 inversion    
    obs_cc = -40.39 
    obs_cc_err = np.array([10, 10]) # obs_cc_err = np.array([-0.06, 4.62])   
    obs_cc_err = obs_cc_err.reshape((2, 1))
    
    # my 11s2 inversion
    obs_m2 = 39.06
    obs_m2_err = np.array([10, 10]) # obs_m2_err = np.array([-7.25E-5, 1.43])
    obs_m2_err = obs_m2_err.reshape((2, 1))
    
    ax.set_ylim([-158.94445, 199.22545])
    if xplot == "R":
        ax.set_xlim([np.min(m1[0]), 0])
    else:
        ax.set_xlim([0, np.max(m1[0])])
    ax.set_ylabel("$c_{00}$ ($\mu$Hz)", fontsize=8, ha="right")
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    
elif yplot == "Q":
    y1  = m1[2]
    y2  = m2[2]
    
    # my 10s2 inversion
    obs_m1 = 670  
    obs_m1_err = None
    
    # my 11s2 inversion
    obs_m2 = 95
    obs_m2_err = None
    
    ax.set_ylabel("$Q$", fontsize=8)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    
elif yplot == "freq":
    y1  = m1[3]
    y2  = m2[3]

    # my 10s2 inversion
    obs_m1 = 4040.09  
    obs_m1_err = None
    
    # my 11s2 inversion
    obs_m2 = 4073.11
    obs_m2_err = None
    
    ax.set_ylabel("$f_c$ ($\mu$Hz)", fontsize=8, va="bottom")
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    
# c00
ax.plot(m1[0], y1, label='${}_{10}S_2$', color='r', 
        linestyle='-', linewidth=0.75)  
ax.errorbar(x=param_m1, y=obs_m1, yerr=obs_m1_err, 
            mew=0.75, color='r',  mfc='w', elinewidth=0.75, mec='r',
            marker='D', markersize=4, capsize=2.5, zorder=50,) 
            
ax.plot(m2[0], y2, label='${}_{11}S_2$', color='k', 
        linestyle='-', linewidth=0.75)
ax.errorbar(x=param_m2, y=obs_m2, yerr=obs_m2_err, 
            mew=1, color='k',  mfc='w', elinewidth=0.75, mec='k',
            marker='D', markersize=4, capsize=2.5, zorder=50,) 
    
if yplot == "c00":
    ax.plot(cc[0], cc[1], label='${}_{10}S_2$-${}_{11}S_2$', color='grey', 
            linestyle='-', linewidth=0.75)
    ax.errorbar(x=param_cc, y=obs_cc, yerr=obs_cc_err, 
                mew=1, color='grey',  mfc='w', elinewidth=0.75, mec='grey',
                marker='D', markersize=4, 
                capsize=2.5, zorder=50,) 

if xplot == "vs":
    ax.axvspan(0.35, 0.75, alpha=0.60, color=color)
    ax.set_xlabel("$\delta v_s/ v_s$ ($\%$)", fontsize=8)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
elif xplot == "vp":
    ax.axvspan(2.5, 4.5, alpha=0.60, color=color)
    ax.set_xlabel("$\delta v_p/ v_p$ ($\%$)", fontsize=8)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
elif xplot == "rho":
    ax.axvspan(18, 37, alpha=0.60, color=color)
    ax.set_xlabel("$\delta \\rho / \\rho$ ($\%$)", fontsize=8)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
elif xplot == "R":
    ax.axvspan(-0.75, -0.35, alpha=0.60, color=color)
    ax.set_xlabel("$\delta R/ R$ ($\%$)", fontsize=8)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))

if yplot == "c00" and xplot == "vs":
    plt.legend(loc="center right", handlelength=0.65, 
                handletextpad=0.15, labelspacing=0.1,
                fontsize=8, frameon=False, bbox_to_anchor=(1.05, 0.42))

elif yplot == "c00" and xplot == "vp":
    plt.legend(loc="upper left", handlelength=0.65, 
                handletextpad=0.15, labelspacing=0.1,
                fontsize=8, frameon=False, bbox_to_anchor=(0, 1.05))
    

# ax.yaxis.set_label_coords(-0.2,0.45)
# ax.xaxis.set_label_coords(0.55,-0.17)

plt.savefig('%s_v_%s.png'%(yplot, xplot),
              dpi=400,
              bbox_inches='tight', 
              pad_inches=0.0,
              )