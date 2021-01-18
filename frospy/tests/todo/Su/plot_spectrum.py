#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:05:03 2019

@author: talavera
"""
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker

#mpl.rcParams['axes.linewidth'] = 1 #set the value globally
##
##event="031111B"
##station="PPTF" #
#event="052413A"
##station="PAYG" #
##event="100494B"
##event="022710A"
##event="060994A"
#station="KBS" #
#
#if station is "PAYG":
#    fw=[4.022, 4.065] # segment fw
#if station is "PPTF":
#    fw=[4.022, 4.055]
#if station is "KBS":
#    fw=[4.022, 4.065]
#    
#home = "/net/home/talavera/eejit/splitting"
#cwd1 = "%s/10s02-11s02/allevents/cst+d00+cc/c00"%home
#cwd2 = "%s/04s00/synthetics/IC/QCC"%home
#s1 = spectrum(data="%s/%s.ahx"%(cwd2,event),
#             syn=[
##                  "%s/%s.PREM.noQCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.s20rts+crust.ahx.syn"%(cwd2,event),
##                  "%s/%s.highVs.noQCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.PREM.QCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.s20rts+crust.QCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.highVs.QCC.ahx.syn"%(cwd2,event),
#                  "%s/%s.mnew.ahx.syn"%(cwd2,event), 
#                  ],
#                  tw=[12,60], fw=fw, # 
#                  minispec=["phase"],# "amp"], 
#                  syn_label=[
##                              "PREM",
##                             "S20RTS+CRUST",
##                             "$2\%$ $V_S$",
#                             "Our new $c_{st}$",
#                             ], 
##                  modes=["10S2", "11S2"],
#                  station=station,
##                  segment_file = "%s/%s/%s.dat"%(cwd1,event,event)
#                  width=2.,
##                  cmap='BlueBlackGreysRed',
#                  )
#
#s2 = spectrum(data="%s/%s.ahx"%(cwd2,event),
#             syn=[
##                  "%s/%s.PREM.noQCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.s20rts+crust.ahx.syn"%(cwd2,event),
##                  "%s/%s.highVs.noQCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.PREM.QCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.s20rts+crust.QCC.ahx.syn"%(cwd2,event),
##                  "%s/%s.highVs.QCC.ahx.syn"%(cwd2,event),
#                  "%s/%s.mnew.ahx.syn"%(cwd2,event), 
#                  ],
#                  tw=[12,60], fw=fw, #
#                  minispec=[ "amp"], 
#                  syn_label=[
##                              "PREM",
##                             "S20RTS+CRUST",
##                             "$+2\%$ $v_s$",
#                             "Our Inversion",
#                             ], 
##                  modes=["10S2", "11S2"],
#                  station=station,
##                  segment_file = "%s/%s/%s.dat"%(cwd1,event,event)
#                  width=2.,
##                  cmap='BlueBlackGreysRed',
#                  )
#
#if event is "052413A":
#    t1 = 'a) 1D Coupling synthetics'
#    t2 = '$M_w$ $8.3$ $\mathrm{Okhotsk,}$ $2013$' 
#    t3 = '$%s$ $\mathrm{station,}$ $12$-$60$$\mathrm{h}$'%station
#    
#
#if event is "031111B":
#    t1 = 'b) 1D Coupling synthetics'
#    t2 = '$M_w$ $9.1$ $\mathrm{Tohoku,}$ $2011$'
#    t3 = '$%s$ $\mathrm{station,}$ $12$-$60$$\mathrm{h}$'%station
#    
#
#if event is "100494B":
#    t1 = 'b) 1D Coupling synthetics'
#    t2 = '$M_w$ $8.4$ $\mathrm{Kuril,}$ $1994$' 
#    t3 = '$%s$ $\mathrm{station,}$ $12$-$60$$\mathrm{h}$'%station
#    
#if event is "022710A":
#    t1 = 'b) 1D Coupling synthetics'
#    t2 = '$M_w$ $8.8$ $\mathrm{Chile,}$ $2010$' 
#    t3 = '$%s$ $\mathrm{station,}$ $12$-$60$$\mathrm{h}$'%station
#    
#
#s1[0].set_size_inches(2.3,0.75)
#s1[0].suptitle('%s\n%s\n%s' % (t1,t2,t3), y=1.6, x=0.53,
#               fontsize=10, weight="bold")
#s1[0].axes[0].tick_params(axis='both', which='both', bottom=False, labelsize=8)
#s1[0].axes[0].spines['bottom'].set_visible(False)
#
#s2[0].axes[0].get_legend().remove()
#s2[0].set_size_inches(2.3,2.25)
#s2[0].axes[0].xaxis.set_major_locator(plt.MaxNLocator(3))
#s2[0].axes[0].tick_params(axis='both', which='major',labelsize=8)
#s2[0].axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
#s2[0].axes[0].set_xlabel("f(mHz)", fontsize=8)
#
#if station is "PAYG":
#    s1[0].axes[0].set_ylabel("Phase", fontsize=8)
#    s2[0].axes[0].set_ylabel("Amplitude", fontsize=8)
##    s2[0].text(0.29,0.75,'${}_{%s}S_{%s}$'%(10,2),fontsize=10)
##    s2[0].text(0.58,0.63,'${}_{%s}S_{%s}$'%(11,2),fontsize=10)
#    s2[0].legend(
#    #            bbox_to_anchor=(0.065, 0.08),
#    #            bbox_to_anchor=(0.05, 0.08),
#    #            loc="upper left", columnspacing=0.75,
#                loc="upper left", bbox_to_anchor=(0.58, 0.98), 
#                frameon=False, borderpad=0,#ncol=5,
#                handlelength=1, handletextpad=0.1,
#                labelspacing=0.2, fontsize=8
#                )
#if station is "PPTF":  
#    s1[0].axes[0].set_yticklabels([])
#    s1[0].axes[0].set_ylabel("", fontsize=8)
#    s2[0].axes[0].set_ylabel("", fontsize=8)
#    s2[0].axes[0].set_yticklabels([])
##    s2[0].text(0.35,0.65,'${}_{%s}S_{%s}$'%(10,2),fontsize=10)
##    s2[0].text(0.74,0.37,'${}_{%s}S_{%s}$'%(11,2),fontsize=10)
#    
#if station in ["SPA","KBS","KEV","ALE"]:  
#    s1[0].axes[0].set_yticklabels([])
#    s1[0].axes[0].set_ylabel("", fontsize=8)
#    s2[0].axes[0].set_ylabel("", fontsize=8)
#    s2[0].axes[0].set_yticklabels([])
#
#s1[0]
#s2[0]
#s1[0].savefig('%s_%s_phase.png'%(event,station),
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.,
#              )
#
#s2[0].savefig('%s_%s_amp.png'%(event,station),
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.,
#              )

#
###------------------------------------------------------------------
#import sys
#import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from matplotlib import ticker
#
#mpl.rcParams['axes.linewidth'] = 1.6 #set the value globally
#
#event="060994A"
#fw=[4.022, 4.06]    
#home = "/net/home/talavera/eejit/splitting"
#cwd1 = "%s/10s02-11s02/allevents/cst+d00+cc/c00"%home
#cwd2 = "%s/04s00/synthetics/IC/QCC"%home
#s3 = spectrum(
##              data="%s/%s.ahx"%(cwd2,event),
#              data="/net/home/talavera/eejit/data/VHZ/%s.ahx"%event,
#              fw=fw, # 
#              seg_reference=True,
#              minispec="amp", #"phase"],# 
##              modes=["10S2", "11S2"],
#              segment_file = "%s/%s/%s.dat"%(cwd1,event,event),
#              width=1.,
##              cmap='GreysRed',
#              runall=True, 
#              plot=True,
#             )
#
#t1 = 'c) Real Spectra data'
#
#if event is "060994A":
#    t2 = '$M_w$ $8.3$  $\mathrm{Bolivia,}$ $1994$' 
#    t3 = ' $\mathrm{All}$ $\mathrm{stations,}$ $12$-$60$$\mathrm{h}$'
#else:
#    t2 = '$\mathrm{%s}$'%event 
#    t3 = ' $\mathrm{All}$ $\mathrm{stations,}$ $12$-$60$$\mathrm{h}$'
#
#s3[0].set_size_inches(2.3,3.32)
#s3[0].suptitle('%s\n%s\n%s' % (t1,t2,t3), y=1.026, x=0.5,
#               fontsize=10, weight="bold")
#s3[0].axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
#s3[0].axes[0].xaxis.set_major_locator(plt.MaxNLocator(3))
#s3[0].axes[0].set_xlabel("f(mHz)", fontsize=8)
#s3[0].axes[0].set_ylabel("", fontsize=8)
#s3[0].axes[0].tick_params(axis='both', which='both', labelsize=8)
#s3[0].axes[0].set_yticklabels([])
#s3[0]
#s3[0].savefig('%s_all_amp.png'%event,
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.,
#              )

##------------------------------------------------------------------
#import numpy as np
#import matplotlib.pyplot as plt
#import sys
#
#event = "031111B"
##cwd = "/net/home/talavera/eejit/splitting/02s12/synthetics/R"
#cwd = "/net/home/talavera/eejit/data/R"
#
#sta= "COLA"
#plots = ["amp"]#, "gcp"]
#ax_dum = []
#
#for plot in plots:
#    dum = spectrum(
#                  #data='%s/%s.ahx'% (cwd, event),
#                  data='%s/%s/%s.ahx'% (cwd, "zero", event),
#                  syn=[
#                       #'%s/%s.ahx.R-0.2.syn' % (cwd, event),
#                       '%s/%s/%s.ahx.syn.fil' % (cwd, "anticorr", event),
#                       '%s/%s/%s.ahx.syn.fil' % (cwd, "corr", event),
#                       '%s/%s/%s.ahx.syn.fil' % (cwd, "zero", event),
#                       ],
#                  syn_label=[
#                             #'R=-0.2',
#                             'Anticorrelated 3D Q',
#                             'Correlated 3D Q',
#                             'Zero 3D Q' ,
#                             ], 
#                  #segment_file='../segments/%s.dat' %event,
#    #              fw=[2.720,2.760], tw=[3,50],
#                  fw=[2.39,2.60], tw=[5,60],
#                  minispec = [plot],
#                  modes=True,
#                  #nowindow=True,
#                  station = sta,
#                  )
#    ax_dum.extend(dum)
#
##------------------------------------------------------------------
#import numpy as np
#import matplotlib.pyplot as plt
#import sys
#event="052413A"
#station="KURK"
#home = "/net/home/talavera/eejit/splitting"
#cwd1 = "%s/05s00/deepevents/self/c00=10"%home
#cwd2 = "%s/05s00-13s02/deepevents/cst+d00+cc/cst_05s00_13s02"%home
#s1 = spectrum(data="%s/%s.ahx"%(cwd1,event),
#             syn=["%s/%s.ahx.syn"%(cwd1,event), 
#                  "%s/%s.ahx.syn"%(cwd2,event)], 
#                  tw=[35,90], fw=[4.820000, 4.901123], 
#                  minispec=["phase"], 
#                  syn_label=["SC, misfit = 0.119", 
#                             "GC, misfit = 0.036"], 
##                  modes=["13S2", "5S0"],#, "1T25", "3T15", "9S7"],
#                  station=station,
#                  segment_file = "%s/%s/%s.dat"%(cwd1,event,event)
#                  )
#
#s = spectrum(data="%s/%s.ahx"%(cwd1,event),
#             syn=["%s/%s.ahx.syn"%(cwd1,event), 
#                  "%s/%s.ahx.syn"%(cwd2,event)], 
#                  tw=[35,90], fw=[4.820000, 4.901123], 
#                  minispec=["amp"], 
#                  syn_label=["SC, misfit = 0.119", 
#                             "GC, misfit = 0.036"], 
#                  modes=["13S2", "5S0"],#, "1T25", "3T15", "9S7"],
#                  station=station,
#                  segment_file = "%s/%s/%s.dat"%(cwd1,event,event)
#                  )
#s[0].axes[0].get_legend().remove()
#s[0].set_size_inches(6,3)
#s[0].legend(#bbox_to_anchor=(0.15, 1.1),
#            bbox_to_anchor=(0.15, 0.08),
#            loc="upper left", columnspacing=0.75,
#            frameon=False, ncol=3,
#            handlelength=1, handletextpad=0.25)
#s[0]

#s[0].savefig('%s_%s_amp.png'%(event,station),
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
##              pad_inches=0.025,
#              )

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker

mpl.rcParams['axes.linewidth'] = 1.8 #set the value globally
#
#event="031111Q"

# 2s0-7s2
#station="RER" #031111Q:CTAO, 060994Q:TRQA, 052413A: RER
#event="052413A"

# 5s0
#station="KEV" #
#event="060994A"

##11s0
station="INCN" #
event="053015A"

if station is "KEV":
#    fw=[4.820000, 4.901123] # segment fw
    fw=[4.873657,4.901123] # segment fw
    tw=[35,90]
    syn_label=[
#               r"PREM",
               r"SC",# $m^{c_{st}}$ = 0.07", 
               r"GC", #$m^{c_{st}}$ = 0.01", 
              ]
    home = "/net/home/talavera/eejit/splitting"
    cwd0 = "%s/05s00/synthetics/PREM"%home
    cwd1 = "%s/05s00/deepevents/self/c00=10"%home
    cwd2 = "%s/05s00-13s02/deepevents/cst+d00+cc/cst_05s00_13s02"%home
    syn=[
#         "%s/%s.ahx.syn"%(cwd0,event), 
         "%s/%s.ahx.syn"%(cwd1,event), 
         "%s/%s.ahx.syn"%(cwd2,event),  
        ]

if station is "RER":
    fw=[2.50, 2.53]
    tw=[5,150]
    syn_label=[
#               r"PREM",
               r"SC",# $m^{c_{st}}$ = 0.33", 
               r"GC", #$m^{c_{st}}$ = 0.01", 
              ]
    home = "/net/home/talavera/eejit/splitting"
    cwd0 = "%s/02s00-07s02/synthetics/PREM"%home
    cwd1 = "%s/02s00-07s02/allevents/c00+d00/c00=-10"%home
    cwd2 = "%s/02s00-07s02/allevents/cst+d00+cc/cst_02s00_c20=15"%home
    syn=[
#         "%s/%s.ahx.syn"%(cwd0,event), 
         "%s/%s.ahx.syn"%(cwd1,event), 
         "%s/%s.ahx.syn"%(cwd2,event),  
        ]

if station is "INCN":
    fw=[9.845,9.920]
    tw=[5,50]
    syn_label=[
#               r"PREM",
               r"GC", #$m^{c_{st}}$ = 0.01", 
              ]
    home = "/net/home/talavera/eejit/splitting"
    cwd0 = "%s/11s00-27s02/synthetics/PREM"%home
    cwd1 = "%s/11s00-27s02/allevents/cst+d00+cc/cst_27s02+c00=20"%home
    syn = [
#           "%s/%s.ahx.syn"%(cwd0,event), 
           "%s/%s.ahx.syn"%(cwd1,event)
           ]
    
s1 = spectrum(data="%s/%s.ahx"%(cwd1,event),
             syn=syn,
                  tw=tw, fw=fw, # 
                  minispec=["phase"],# "amp"], 
                  syn_label=syn_label, 
#                  modes=["13S2", "5S0"],#, "1T25", "3T15", "9S7"],
                  station=station,
#                  segment_file = "%s/%s/%s.dat"%(cwd1,event,event)
                  width=1.,
                  cmap='BlackGreysRed', #'Grays',
                  )

s2 = spectrum(data="%s/%s.ahx"%(cwd1,event),
             syn=syn,
                  tw=tw, fw=fw, # 
                  minispec=["amp"], 
#                  minispec="amp", 
                  syn_label=syn_label,
#                  modes=["13S2", "5S0"],#, "1T25", "3T15", "9S7"],
                  station=station,
#                  segment_file = "%s/%s/%s.dat"%(cwd2,event,event),
                  width=1.,
                  cmap='BlackGreysRed', #'Grays',
                  )

if station is 'RER':
#    t1 = 'Radial mode ${}_{%s}S_{%s}$' % (2, 0)
#    t2 = '$M_w$ $9.1$ $\mathrm{Tohoku,}$ $2011$'
    t2 = '$M_w$ $8.3$ $\mathrm{Okhotsk,}$ $2013$' 
    t3 = '$%s$ $\mathrm{station,}$ $%s$-$%s\mathrm{h}$'%(station, tw[0], tw[1])
#    t4 = 'a) PREM synthetics: ${}_{%s}S_{%s}$' % (2, 0)
    t4 = 'd) Our inversion: ${}_{%s}S_{%s}$' % (2, 0)
    
if station is 'KEV':
#    t1 = 'Radial mode ${}_{%s}S_{%s}$' % (5, 0)
    t2 = '$M_w$ $8.3$  $\mathrm{Bolivia,}$ $1994$' 
    t3 = '$%s$ $\mathrm{station,}$ $%s$-$%s\mathrm{h}$'%(station, tw[0], tw[1])
#    t4 = 'b) PREM synthetics: ${}_{%s}S_{%s}$' % (5, 0)
    t4 = 'e) Our inversion: ${}_{%s}S_{%s}$' % (5, 0)
    
if station is 'INCN':
#    t1 = 'Radial mode ${}_{%s}S_{%s}$' % (11, 0)
    t2 = '$M_w$ $7.9$ $\mathrm{Bonin}$ $\mathrm{Islands,}$ $2015$' 
    t3 = '$%s$ $\mathrm{station,}$ $%s$-$%s\mathrm{h}$'%(station, tw[0], tw[1])
#    t4 = 'c) PREM synthetics: ${}_{%s}S_{%s}$' % (11, 0)
    t4 = 'f) Our inversion: ${}_{%s}S_{%s}$' % (11, 0)

s1[0].set_size_inches(2.3,0.75)
s1[0].suptitle('%s\n%s\n%s' % (t2,t3,t4), y=1.6, x=0.53,
               fontsize=10, weight="bold")
s1[0].axes[0].tick_params(axis='both', which='both', bottom=False, labelsize=8)
s1[0].axes[0].spines['bottom'].set_visible(False)

s2[0].axes[0].get_legend().remove()
s2[0].set_size_inches(2.3,2.25)
s2[0].axes[0].xaxis.set_major_locator(plt.MaxNLocator(3))
s2[0].axes[0].tick_params(axis='both', which='major',labelsize=8)
s2[0].axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
s2[0].axes[0].set_xlabel("f(mHz)", fontsize=8)

if station is "RER":
    s1[0].axes[0].set_ylabel("Phase", fontsize=8)
    s2[0].axes[0].set_ylabel("Amplitude", fontsize=8)
    s2[0].text(0.24,0.60,'${}_{%s}S_{%s}$'%(2,0),fontsize=10)
    s2[0].text(0.55,0.35,'${}_{%s}S_{%s}$'%(7,2),fontsize=10)
    s2[0].legend(
    #            bbox_to_anchor=(0.065, 0.08),
    #            bbox_to_anchor=(0.05, 0.08),
    #            loc="upper left", columnspacing=0.75,
                loc="upper left", bbox_to_anchor=(0.55, 0.95), 
                frameon=False, borderpad=0,#ncol=5,
                handlelength=1, handletextpad=0.1,
                labelspacing=0.2, 
                fontsize=8
                )
if station is "KEV":  
    s1[0].axes[0].set_yticklabels([])
    s1[0].axes[0].set_ylabel("", fontsize=8)
    s2[0].axes[0].set_ylabel("", fontsize=8)
    s2[0].axes[0].set_yticklabels([])
    s2[0].text(0.43,0.75,'${}_{%s}S_{%s}$'%(5,0),fontsize=10)
#    s2[0].text(0.74,0.37,'${}_{%s}S_{%s}$'%(7,2),fontsize=10)
if station is "INCN":  
    s1[0].axes[0].set_yticklabels([])
    s1[0].axes[0].set_ylabel("", fontsize=8)
    s2[0].axes[0].set_ylabel("", fontsize=8)
    s2[0].axes[0].set_yticklabels([])
    s2[0].text(0.76,0.37,'${}_{%s}S_{%s}$'%(11,0),fontsize=10)
    s2[0].text(0.39,0.75,'${}_{%s}S_{%s}$'%(27,2),fontsize=10)
#    s2[0].text(0.74,0.37,'${}_{%s}S_{%s}$'%(7,2),fontsize=10)
s1[0]
s2[0]
#s1[0].savefig('%s_%s_phase_measured.png'%(event,station),
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.,
#              )
#
#s2[0].savefig('%s_%s_amp_measured.png'%(event,station),
##              orientation='landscape', 
#              dpi=400,
#              bbox_inches='tight', 
#              pad_inches=0.,
#              )