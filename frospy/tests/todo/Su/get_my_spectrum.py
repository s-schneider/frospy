#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:26:29 2018

@author: talavera
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from frospy.util.read import read as read_modes
from frospy.spectrum.controllers import (read_data, read_syn,
                                       _plot_modes_in_ax, get_mode)
from frospy.util.base import (fourier_transform, get_times, taper)
from frospy.core.segment import read as read_segments

#event='060994A'
event='031111B'
data='PREM/%s.ahx'%event
syn=[
     'PREM/%s.ahx.syn'%event, # self
#     'PREM/%s.ahx.prem.syn'%event, # coupled 
#     'IC/%s.zero.ahx.syn.fil'%event, # group S20RTS
#     'S20RTS/%s.ahx.s20rts.syn'%event,
#     'IC/%s.woodhouse.ahx.syn.fil'%event,
#     'IC/%s.romanowicz.ahx.syn.fil'%event,
#     'IC/%s.tromp.ahx.syn.fil'%event,
#     'IC/%s.begtramp.ahx.syn.fil'%event,
#     'QM1/%s.ahx.syn'%event,
#     'DE/%s.ahx.syn'%event,
#     'SS/%s.ahx.syn'%event,
     'HT/%s.ahx.syn'%event,
#     'Qtest/%s.ahx.syn'%event,
     'inv-self/%s.ahx.syn'%event,
#     'inv-cross/%s.ahx.syn'%event
     ] 

syn_label=[
           'PREM',
#           'PREM-cross', 
#           'S20RTS-group', 
#           'S20RTS',
#           'IC Woodhouse', 
#           'IC Romanowicz', 
#           'IC Tromp', 
#           'IC BegTramp',
#           'Widmer et al, 1990',
#           'Durek & Ekstrom, 1995', 
#           'Masters & Widmer (REM)', 
           'He & Tromp, 1996', 
#           'Qtest',
           'My Inversion'
#           'My Inversion: self-coupling',
#           'My Inversion: cross-coupling'
           ]

segment = read_segments('segments/%s.dat'%event)
stations = [s.station for s in segment]

#mode=['${}_1 S_0$']; fw = [1.618, 1.641]; tw = [5,280]
#mode=['${}_2 S_0$']; fw = [2.489, 2.534]; tw = [30,150]
#mode=['${}_3 S_0$']; fw = [3.255, 3.285]; tw = [35,105]
mode=['${}_4 S_0$']; fw = [4.085, 4.125]; tw = [25,100]
#mode=['${}_5 S_0$']; fw = [4.872, 4.895]; tw = [35,90]
#mode=['${}_6 S_0$']; fw = [5.724, 5.759]; tw = [20,75]

st, st_work, syn = read_data(data, syn)
st_syn, st_syn_work = read_syn(syn, st, syn_label) 

taper_shape='hanning'
Htstart = tw[0] * 3600.
Htend = tw[1] * 3600.
wstart = fw[0]
wend = fw[1]
Htstart_org, Htend_org = Htstart, Htend
#color = plt.cm.get_cmap("gist_rainbow", len(st_syn))
#color = ['grey', 'magenta', 'limegreen', 'c', 'b', 'r']; 
color = ['grey', 'b', 'r']; 
#color = ['c', 'b']; 
#color = ['grey', 'r']; 

i = 0; j = 0; key = ['nps']
#while key != 'quit':
#for i in range(1):
for i in [21]: #i=station
        # startup
        Fxx = []; Fxxsyn = []
        xp = []; xpsyn = []
        freq = [] 

        # data   
        tr_data = st[i]
        times = get_times(tr_data, Htstart, Htend, Htstart_org, Htend_org)
        tstart, nstart, nend, t, Htstart, Htend, Terr = times
        trfill = taper(nstart, nend, tr_data, taper_shape) # taper
        # Fouriertransform and filter
        f, FT, XP, delomeg = fourier_transform(tr_data, Htstart, Htend, nstart, 
                                               nend, taper_shape, 1) 
        freq.append(f)
        Fxx.append(FT)
        xp.append(XP)
    
        # syn
        for syn in st_syn:
            tr_syn = syn[0]
            trsynfill = taper(nstart, nend, tr_syn, taper_shape)
            void, FT, XP, void = fourier_transform(tr_syn, Htstart, Htend, 
                                                   nstart, nend, taper_shape, 1)
            Fxxsyn.append(FT)
            xpsyn.append(XP)
        
        #setting limits
        startlabel = []; endlabel = [];
        startlabel.append(int(np.round(wstart * 2. * np.pi/(1000. * delomeg))) + 1)
        endlabel.append(int(np.round(wend * 2. * np.pi/(1000. * delomeg))) + 1)
        startlabel=startlabel[0]
        endlabel=endlabel[0]        
        
        #plotting
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(2,1,1)
        sta = '%s' % (tr_data.stats.station)
        ax.plot(f[startlabel:endlabel+1], xp[0][startlabel:endlabel+1],
                linestyle='solid', color='black', linewidth=3)
        for k in range(len(Fxxsyn)):
#        for k in range(1):
            ax.plot(f[startlabel:endlabel+1], xpsyn[k][startlabel:endlabel+1],
                    linestyle='dashed', color=color[k], linewidth=3.5)
        box = ax.get_position() # Shrink current axis's height by 35% on the bottom
        ax.set_position([box.x0, box.y0 + box.height * 0.65,
                             box.width, box.height * 0.35])
        ax.set_ylabel('Phase', fontsize=20)
        ax.get_xaxis().set_ticklabels([])
        ax.set_xlim(f[startlabel], f[endlabel+1])
        ax.set_ylim(-3.2,3.2)
        ax.tick_params(axis = 'both', which = 'major', labelsize = 13, zorder=0)
        ax.set_title(r'Station %s, %s - %s' % (sta, tw[0], tw[1]), fontsize=25)
    
    
        ax = fig.add_subplot(2,1,2)
        ax.plot(f[startlabel:endlabel+1], abs(Fxx[0][startlabel:endlabel+1]),
                linestyle='solid', color='black', linewidth=3, label='data')
        for k in range(len(Fxxsyn)):
#        for k in range(1):
            ax.plot(f[startlabel:endlabel+1], abs(Fxxsyn[k][startlabel:endlabel+1]),
                    linestyle='dashed', color=color[k], 
                    linewidth=3.5, label=syn_label[k])
        box = ax.get_position() # Shrink current axis's height by 10% on the bottom
        ax.set_position([box.x0, box.y0 + box.height * 0.15,
                         box.width, box.height * 1.5])
        modedata = read_modes()#; modes = modedata
        #modes=get_mode(['3s0', '9S2', '8S2'], modedata)
        f_lim = [f[startlabel], f[endlabel+1]]
        _plot_modes_in_ax(ax, modedata, f_lim)
        ax.set_ylabel('Amplitude', fontsize=20)
        ax.set_xlabel('frequency (mHz)', fontsize=20)
#        ax.set_ylim(0,5.1e-7)
        ax.set_ylim(0,1e-6)
        ax.tick_params(axis = 'both', which = 'major', labelsize = 13, zorder=0)
        ax.legend(fontsize=16, frameon=False, loc='upper left')
        mf = ticker.ScalarFormatter(useMathText=True)
        mf.set_powerlimits((-2,2)) # for exponents greater than +-2
        plt.gca().yaxis.set_major_formatter(mf)
        plt.show()
#        plt.savefig('%s-%s-amp'%(event,sta), bbox_inches='tight', 
#                    pad_inches=0, dpi=350)
    
#        key = raw_input('Please an action (nsp, bts, psp, tw, fw, quit): ')
#        if i >= len(st_work): i = 0
#        elif key == 'next' or key == 'nsp' or key == "": i += 1
#        elif key == 'quit' or key == 'exit' or key == 'q': break
#        elif key == 'bts': i = 0
#        elif key == 'fw':
#            fw1 = raw_input('Input start fw: '); wstart=float(fw1)
#            fw2= raw_input('Input end fw: '); wend=float(fw2)
#        elif key == 'tw':
#            tw1 = raw_input('Input start tw: '); tw[0]=float(tw1)
#            tw2 = raw_input('Input end tw: '); tw[1]=float(tw2)
#            Htstart = tw[0] * 3600.
#            Htend = tw[1] * 3600.
#            Htmid = (tw[1] - tw[0])/2.
#        elif key == 'psp':
#            if i != 0:
#                i -= 1
#            else:
#                i = 0
#                print('Already at the beginning.\n')
