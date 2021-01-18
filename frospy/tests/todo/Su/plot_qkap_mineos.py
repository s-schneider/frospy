#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:04:50 2019

@author: talavera
"""
import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, basename, splitext
from frospy.util.base import sort_human

greens = plt.cm.get_cmap("Greens", 50)
blues = plt.cm.get_cmap("Blues", 50)

#colors = [
##         blues(10), blues(7), blues(5), #blues(1), 
##          blues(15), blues(12), blues(9), blues(6), blues(3), blues(1), 
#          blues(50), blues(45), blues(40), blues(35), blues(30), blues(25), blues(20), 
##         greens(2), greens(4), greens(6), greens(8), greens(10),greens(12),]
##         greens(8), greens(10),greens(12),]
##         greens(2), greens(4), greens(6), greens(8), 
#        ]

Q_PREM =  [5327,1499,1242,1083,969,921,913,881,852,840,829,832]
Q_CC =  [5983,1856,1788,1228,1181,1006,1083,905,1028,937,0, 985]
err_u = [40, 144, 113, 318,  85,  19,  103, 121,75,  84, 0, 58]
err_l = [127, 143, 122, 35,  34,  8,   43,  7,  33,  45, 0, 49]
data = []
for Q, u, l in zip(Q_CC, err_u, err_l):
    eval = np.array([l, u]).reshape((2, 1))
    data.append([Q, eval])

Q_QM1 = [5724.40,1848.47,1801.83,1410.45,1215, 1250,  1111.11,0,1111, 1219,0,0]
err_u = [-262,   225.35, 485.42, 142.68, 60.69,126.26,161,    0,226.14,316,0,0]
err_l = [-262,   225.35, 485.42, 142.68, 60.69,126.26,161,    0,226.14,316,0,0]
QM1 = []
for Q, u, l in zip(Q_QM1, err_u, err_l):
    eval = np.array([l, u]).reshape((2, 1))
    QM1.append([Q, eval])

   
dhome = "/net/home/talavera/codes/mineos/DEMO/Qkappa/"
mdir = [
        "prem_ocean_S_ICQk*",
        "prem_ocean_S_OCQk*",
        "prem_ocean_S_LMQk*",
        "prem_ocean_S_UMQk*",
        "prem_ocean_S_ATQk*",
        "prem_ocean_S_LTQk*",
        "prem_ocean_S_UALQk*",
        "prem_ocean_S_UMM*",
#        "prem_ocean_S_allQk*",
        ]

#yaxs_split=None
                        
ifiles = [glob.glob(join(dhome,d)) for d in mdir]

for A in ifiles: 
    A = sort_human(A)
    n = A[0].split("/")[-1].split("_")[-2]
    # first color is black and assign later
    colors = [None,blues(45),blues(40),blues(35),blues(30),blues(25),]

    if n == "ICQk":
        for i,a in enumerate(A): # best value alwas in zero position
            if a.find("1200") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-100,800,1700,2300,]
        gridspec_kw = {'height_ratios': [0.3,1]}
#        yaxs_split=None
    elif n == "OCQk":
        for i,a in enumerate(A):
            if a.find("18000") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-4500,-1000,-600,2000,]
        gridspec_kw = {'height_ratios': [1,0.3]}
#        yaxs_split=None
    elif n == "LMQk":
        for i,a in enumerate(A):
            if a.find("18000") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-4500,-1000,-250,2000,]
        gridspec_kw = {'height_ratios': [1,0.3]}
#        yaxs_split=None
    elif n == "UMQk":
        for i,a in enumerate(A):
            if a.find("1100") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-50,800,1750,2250,]
        gridspec_kw = {'height_ratios': [0.3,1]}
#        yaxs_split=None
    elif n == "ATQk":
        for i,a in enumerate(A):
            if a.find("140") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-150,700,1900,2300,]
        gridspec_kw = {'height_ratios': [0.3,1]}
        colors = [None,blues(50),blues(45),blues(40),blues(35),blues(30),blues(25)]
#        yaxs_split=None
    elif n == "LTQk":
        yaxs_split=[-50,750,1300,2300,]
        gridspec_kw = {'height_ratios': [0.3,1]}
        colors = [None,blues(50),blues(45),blues(40),blues(35),blues(30),blues(25)]
#        yaxs_split=None
    elif n == "UALQk":
        for i,a in enumerate(A):
            if a.find("1300") > 0: 
                A.insert(0, A.pop(i))
        yaxs_split=[-50,700,1700,2250,]
        gridspec_kw = {'height_ratios': [0.3,1]}
#        yaxs_split=None
    elif n == "allQk":
        yaxs_split=None
    else:
        colors = ["b", "k", "grey","lightgrey"]
        yaxs_split=[-1700,-1500,-50,800,]
        gridspec_kw = {'height_ratios': [1,0.3]}
        
    #splitaxis
    if isinstance(yaxs_split, list):
        fig, (axt, ax) = plt.subplots(2,1,sharex=True,
                                        gridspec_kw=gridspec_kw,
                                        facecolor='w')
        ax.set_ylim(yaxs_split[0], yaxs_split[1])
        axt.set_ylim(yaxs_split[2], yaxs_split[3])
        
        ax.spines['top'].set_visible(False)
        axt.spines['bottom'].set_visible(False)
        ax.tick_params(top=False)
        axt.tick_params(bottom=False)

        # directly from:
        # https://matplotlib.org/examples/pylab_examples/broken_axis.html
        d = .01  # how big to make the diagonal lines in ax
        kwargs = dict(transform=ax.transAxes, color='k',
                                            clip_on=False)
        axt.plot((-d,+d),(1-d,1+d),**kwargs) # bottom/left
        axt.plot((1-d,1+d),(1-d,1+d),**kwargs) # bottom/right
        
        kwargs.update(transform=axt.transAxes, color='k',
                                              clip_on=False)
        ax.plot((1-d,1+d),(-d,+d),**kwargs) # top/right
        ax.plot((-d,+d),(-d,+d),**kwargs) # top/left
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    fig.set_size_inches(6, 4)

    for j, f in enumerate(A):
        Qk = float(f.split("/")[-1].split("_")[-1])
        
        try:
            m = f.split("/")[-1].split("_")[-2].split("UMM")[1]
        except Exception:
            m = f.split("/")[-1].split("_")[-2].split("Qk")[0]
            mq = f.split("/")[-1].split("_")[-1]
        
#        if m == "QL6":
#            name = "%s $Q_\kappa^{UM}=$%s" % (m, int(Qk))
        if m in ["LT", "AT", "UM", "LM", "OC", "IC", "UAL"]:
            name = "$Q_\kappa=$%s" % (mq)
        else:
            name = "%s" % (m)
        
        f = open(f, 'r') 
        lines = f.read().splitlines()
        lines = lines[-12:]
        color = colors[j]
            
        for i, (l, Q0, d, q) in enumerate(zip(lines, Q_PREM, data, Q_QM1)):
#       for i, (l, Q0, d) in enumerate(zip(lines, Q_PREM, data)):
            l = l.split()
            mode = '%s%s%s' % (l[0], l[1], l[2])
            Qs = int(float(l[-2]))
            if i == 0 and j == 1:
                #measured
                ax.errorbar(i, d[0]-Q0, yerr=d[1], label="GC",
                            linestyle='None',color="r", ms=7, 
                            marker='o', capsize=3,zorder=50)
                if isinstance(yaxs_split, list):
                    axt.errorbar(i,d[0]-Q0,yerr=d[1],ms=7,linestyle='None', 
                                 color="r", marker='o', capsize=3,zorder=50)
#               #QM1
                ax.scatter(i, q-Q0, label="QM1$^{obs}$", 
                           color="gray", marker='D')
                if isinstance(yaxs_split, list):
                    axt.scatter(i,q-Q0,color="gray",marker='D')
            if i == 0: 
                if ((m == "IC" and int(mq) == 1200) or
                    (m == "OC" and int(mq) == 18000) or
                    (m == "LM" and int(mq) == 18000) or
                    (m == "UM" and int(mq) == 1100) or
                    (m == "LT" and int(mq) == 50) or
                    (m == "AT" and int(mq) == 140) or
                    (m == "UAL" and int(mq) == 1300)):
                    color='k'
                    zo = 100
                else:
                    zo = 1
                        
                # synthetic
                ax.scatter(i,Qs-Q0, marker='s',color=color,s=25,
                           zorder=zo,label=name,)
                if isinstance(yaxs_split, list):
                    axt.scatter(i,Qs-Q0,marker='s',color=color,s=25,
                                zorder=zo,label=name)
                    
            if i != 10:
                if i == 11:
                    i = 10
                        
                if ((m == "IC" and int(mq) == 1200) or
                    (m == "OC" and int(mq) == 18000) or
                    (m == "LM" and int(mq) == 18000) or
                    (m == "UM" and int(mq) == 1100) or
                    (m == "LT" and int(mq) == 50) or
                    (m == "AT" and int(mq) == 140) or
                    (m == "UAL" and int(mq) == 1300)):
                    color='k'
                    zo = 100
                else:
                    zo = 1
                        
                # synthetic
                ax.scatter(i,Qs-Q0, marker='s',color=color,s=25,
                           zorder=zo)
                if isinstance(yaxs_split, list):
                    axt.scatter(i,Qs-Q0, marker='s',color=color,s=25,
                                zorder=zo)
                    
                # measured
                ax.errorbar(i, d[0]-Q0, yerr=d[1], linestyle='None',
                            color="r", marker='o', ms=7,capsize=3,zorder=50)
                if isinstance(yaxs_split, list):
                    axt.errorbar(i, d[0]-Q0, yerr=d[1], linestyle='None',
                                color="r", marker='o', ms=7,capsize=3,zorder=50)
                    
                if i != 7 and i != 10 and i != 11: #QM1
                    ax.scatter(i, q-Q0, color="gray", marker='D')
                    if isinstance(yaxs_split, list):
                        axt.scatter(i, q-Q0, color="gray", marker='D')
   
#    modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,12)]
    modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
    modename.append('${}_{11}S_0$')

    if m not in ["LT", "AT", "UM", "LM", "OC", "IC", "UAL"]:
        m = "h) 1D Q models"
    elif m == "IC":
        m = "a) Inner Core"
        ax.set_ylabel("$\delta Q$")
    elif m == "OC":
        m = "b) Outer Core"
    elif m == "LM":
        m = "c) Lower Mantle"
        axt.set_ylabel("$\delta Q$")
    elif m == "UM":
        m = "d) Upper Mantle (220-660 km)"
    elif m == "AT":
        m = "e) Lower Lithosphere (80-220 km)"
        ax.set_ylabel("$\delta Q$")
    elif m == "LT":
        m = "f) Upper Lithosphere (24-80 km)"
    elif m == "UAL":
        m = "g) Upper Mantle and Lithosphere (24-660 km)"
        ax.set_ylabel("$\delta Q$")
        
    if isinstance(yaxs_split, list):
        axt.axhline(y=0, color="silver", zorder=0)
        ax.axhline(y=0, color="silver", zorder=0)
        axt.set_title("%s" % m, weight="bold")
    else:
        ax.axhline(y=0, color="silver", zorder=0)
        ax.set_title("%s" % m, weight="bold")
    #    ax.set_ylim(-10,800)
    
    ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
    ax.set_xticklabels(modename, rotation='horizontal', fontsize=12)
#    ax.set_xlim(-0.5,11.5)
    ax.set_xlim(-0.5,10.5)
    # removing whiskers from legend
    ehandles0, elabels = ax.get_legend_handles_labels()
    ehandles = [h for h in ehandles0[0:-1]]
    ehandles.append(ehandles0[-1][0])
    
    # QM1 exception
    try:
        ehandles.insert(0,ehandles.pop(elabels.index("QM1$^{obs}$")))
        elabels.insert(0,elabels.pop(elabels.index("QM1$^{obs}$")))
#        ehandles.insert(1,ehandles.pop(elabels.index("QM1")))
#        elabels.insert(1,elabels.pop(elabels.index("QM1")))
    except Exception:
        pass

    ehandles.insert(0,ehandles.pop(elabels.index("GC")))
    elabels.insert(0,elabels.pop(elabels.index("GC")))  
    
    if yaxs_split is None:
#        ax.set_ylim(-50,800)
        ax.legend(ehandles, elabels, ncol=3, loc="upper right",
                  handletextpad=0.2, columnspacing=0.3,
                  frameon=False, borderaxespad=0.)
    else:
        if n in ["LTQk", "ATQk", "UMQk", "OCQk", "ICQk", "UALQk", "LMQk"]:
            axt.legend(ehandles, elabels, ncol=3, loc="upper right",
                      handletextpad=0.2, columnspacing=0.3,
                      frameon=False, borderaxespad=0.)
        else:
            ax.legend(ehandles, elabels, ncol=3, loc="lower right",
                      handletextpad=0.2, columnspacing=0.3,
                      frameon=False, borderaxespad=0.)
    plt.show()        
#    fig.savefig('Qk_%s.png'%n,orientation='landscape', 
#                dpi=400, bbox_inches='tight', pad_inches=0.02)