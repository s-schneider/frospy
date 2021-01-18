#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:26:29 2018

@author: talavera
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
from frospy.util.read import read_modes
from frospy.util.read import read_std_cat
from frospy.preprocessing.spectrum import read_data, read_syn, get_times,taper, fourier_transform
from frospy.preprocessing.spectrum import _plot_modes_in_ax, get_mode
from frospy.preprocessing.data_correction import get_baz_synseis
from frospy.preprocessing.spectrum import read_segment_file
from obspy.geodetics import gps2dist_azimuth
#from operator import itemgetter

def misfit(X, Xsyn, s, e):
    d2 = sum(abs(X[s:e+1])**2.)
    rmf = sum((abs(X[s:e+1]) - abs(Xsyn[s:e+1]))**2.) / d2
    cmf = sum(abs(X[s:e+1] - Xsyn[s:e+1])**2.) / d2

    return rmf, cmf # real and complex

event='060994A'
data='PREM/%s.ahx'%event
syn=['PREM/%s.ahx.syn'%event,
     'PREM/%s.ahx.prem.syn'%event, 
     'IC/%s.zero.ahx.syn.fil'%event, # group S20RTS
     'S20RTS/%s.ahx.s20rts.syn'%event,
     'IC/%s.woodhouse.ahx.syn.fil'%event,
     'IC/%s.romanowicz.ahx.syn.fil'%event,
     'IC/%s.tromp.ahx.syn.fil'%event,
     'IC/%s.begtramp.ahx.syn.fil'%event,
     'Qtest/%s.ahx.syn'%event,
     'HT/%s.ahx.syn'%event,
     'DE/%s.ahx.syn'%event,
     'QM1/%s.ahx.syn'%event] 
syn_label=['PREM-self','PREM-cross', 'S20RTS-group', 'S20RTS',
           'IC Woodhouse', 'IC Romanowicz', 'IC Tromp', 'IC BegTramp',
           'Qtest','HT', 'DE', 'QM1']
segment = read_segment_file('segments/%s.dat'%event)

#mode=['${}_1 S_0$']; fw = [1.618, 1.641]; tw_range = range(10,60,5)
#mode=['${}_2 S_0$']; fw = [2.489, 2.534]; tw_range = range(10,60,5)
mode=['${}_3 S_0$']; fw = [3.255, 3.281]; tw_range = range(10,60,5)
#mode=['${}_4 S_0$']; fw = [4.085, 4.125]; tw_range = range(10,60,5)
#mode=['${}_5 S_0$']; fw = [4.872, 4.895]; tw_range = range(10,60,5)
#mode=['${}_6 S_0$']; fw = [5.724, 5.759]; tw_range = range(10,60,5)

cat = read_std_cat(cmt_id=event)
st, st_work, syn = read_data(data, syn)
st_syn, st_syn_work = read_syn(syn, st, syn_label)
 
tw_shift = []; j = 0   

for t in tw_range:
    tw=[t,110]
    taper_shape='hanning'
    Htstart = tw[0] * 3600.
    Htend = tw[1] * 3600.
    Htmid = (tw[1] - tw[0])/2.    
    Htstart_org, Htend_org, Htmid_org = Htstart, Htend, Htmid

    for sta in segment: 
        # startup
        Fxx = []; Fxxsyn = []
        xp = []; xpsyn = []
        freq = [] 
        rmf = []; cmf = []
        
        # from segment file
        wstart = segment[sta]['fw1']
        wend = segment[sta]['fw2']
        xscale = segment[sta]['weight']

        # data   
        tr_data = st.select(station=sta)[0]
        times = get_times(tr_data, Htstart, Htend, Htstart_org, Htend_org)
        tstart, nstart, nend, t, Htstart, Htend, Htmid, Terr = times
        trfill = taper(nstart, nend, tr_data, taper_shape) # taper
        f, FT, XP, delomeg = fourier_transform(tr_data, Htstart, Htend, nstart, nend, taper_shape, 1) # Fouriertransform and filter
        freq.append(f)
        Fxx.append(FT)
        xp.append(XP)
        
        #setting limits
        startlabel = []; endlabel = [];
        startlabel.append(int(np.round(wstart * 2. * np.pi/(1000. * delomeg))) + 1)
        endlabel.append(int(np.round(wend * 2. * np.pi/(1000. * delomeg))) + 1)
        startlabel=startlabel[0]
        endlabel=endlabel[0]
        
        # syn
        for k in range(len(st_syn)):
            tr_syn = st_syn[k].select(station=sta)[0]
            trsynfill = taper(nstart, nend, tr_syn, taper_shape)
            void, FT, XP, void = fourier_transform(tr_syn, Htstart, Htend, nstart, nend, taper_shape, 1)
            Fxxsyn.append(FT)
            xpsyn.append(XP)
            r, c = misfit(Fxx[0]/xscale, Fxxsyn[k]/xscale, startlabel, endlabel)
            rmf.append((1-r)*100); cmf.append((1-c)*100)
            
        #station metadata
        s_lat = tr_data.stats.ah.station.latitude
        s_lon = tr_data.stats.ah.station.longitude
        e_lat = cat[0].origins[0].latitude
        e_lon = cat[0].origins[0].longitude
        baz = get_baz_synseis(e_lat, e_lon, s_lat, s_lon)
        d_km = gps2dist_azimuth(s_lat, s_lon, e_lat, e_lon)[0]/1000.

        tw_shift.append( [[] for k in range((5+len(st_syn)))] )
        tw_shift[j][0].append(sta) # station
        tw_shift[j][1].append(s_lat); tw_shift[j][1].append(s_lon) # lat, lon
        tw_shift[j][1].append(baz); tw_shift[j][1].append(d_km) # back azimuth and epicentral distance
        tw_shift[j][2].append(tw[0]); tw_shift[j][2].append(tw[1]) # time window
        tw_shift[j][3].append(rmf); tw_shift[j][3].append(cmf) # real and complex misfit
        tw_shift[j][4].append(Fxx[0][startlabel:endlabel+1]) # Data
        for k in range(0,len(st_syn)): tw_shift[j][5+k].append(Fxxsyn[k][startlabel:endlabel+1]) # synthetics
        j = j + 1  
        
tw_shift.sort(key=lambda elem: elem[0]) # sort by station name # sort(key=itemgetter(0))

tw_shift.sort(key=lambda elem: elem[0]) # sort by station name # sort(key=itemgetter(0))
stations= [ts[0][0] for ts in tw_shift[::len(tw_range)]]
color = plt.cm.get_cmap("rainbow", len(st_syn))

#-------- all stations  in Amp vs tw plot
seg_dta = []
for sta in stations: # match stations in segments to stations in data
    seg_dta = seg_dta + [ts for ts in tw_shift if ts[0][0]==sta]
    
j=0; k=len(tw_range) # ploting every station tw's as one line
for sta in stations:
    plt.plot([ts[2][0] for ts in tw_shift][j:k], [np.max(np.abs(ts[4])) for ts in tw_shift][j:k], 
             'ko-', linewidth=0.75)
    plt.text(tw_range[0]-0.5, [np.max(np.abs(ts[4])) for ts in tw_shift][j:k][0], sta)
    for i in range(0,len(st_syn)):
        plt.plot([ts[2][0] for ts in tw_shift][j:k], [np.max(np.abs(ts[i+5])) for ts in tw_shift][j:k], 
                 'o--', color=color(i), linewidth=1)
    j=j+len(tw_range); k=k+len(tw_range)
    #plt.show()

i=0; legend=[]
for label in syn_label:
    legend.append(mpatches.Patch(color=color(i), label=label))
    i=i+1
    
plt.ylabel('Amplitude')
plt.xlabel('tw2 (hrs)')
#plt.xlim(min(tw_shift[0][1])-10, max(tw_shift[0][1])+10)
mf = ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2))
plt.gca().yaxis.set_major_formatter(mf)
plt.title('%s with tw2=%s hrs'%(mode[0], tw_shift[0][2][1]))
plt.legend(handles=legend)
plt.show()

#--------- Amp vs tw, one station
station = 'ALE'
one_sta = [ts for ts in tw_shift if ts[0][0]==station]

fig = plt.figure(figsize=(12,5))
ax = fig.add_subplot(1,1,1)
for ts in one_sta: # modes
    for i in range(0,len(st_syn)): 
        shift=np.argmax(abs(ts[5][0]))-np.argmax(abs(ts[i+5][0])) # only works if syn 5 is PREM-self
        plt.plot([x+ts[2][0]-10+shift for x in range(len(ts[i+5][0]))], abs(ts[i+5][0]), 
                  '--', color=color(i), linewidth=0.5)
    shift=np.argmax(abs(ts[4][0]))-np.argmax(abs(ts[5][0]))
    plt.plot([x+ts[2][0]-10-shift for x in range(len(ts[4][0]))], abs(ts[4][0]), 
              'k-', linewidth=0.75)

plt.plot([ts[2][0] for ts in one_sta], [max(abs(ts[4][0])) for ts in one_sta], 
        'ko-', label='data', linewidth=1.5)
for i in range(0,len(st_syn)): 
    plt.plot([ts[2][0] for ts in one_sta], [max(abs(ts[i+5][0])) for ts in one_sta], 
             'o--', color=color(i), label=syn_label[i], linewidth=1.5)
mf = ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2))
plt.gca().yaxis.set_major_formatter(mf)
plt.ylabel('Amplitude')
plt.xlabel('tw2 (hrs)')
plt.title('%s with $tw_2$=%s hrs in %s'%(mode[0], one_sta[0][2][1], station))
plt.legend()
plt.show()

#----- colormap data/syn ratio
#new_seg_sta = seg_sta[:]
#for rm in ["YSS", "HNR", "SUR", "CMO", "WLF", "HAL"]: 
#    new_seg_sta.remove(rm)
#range_ind = [ts[0][0] for ts in tw_shift].index('KIP')

sorting=0  #sorting=0  by latitude; sorting=2 by baz; sorting=3 by epicentral distance 
tw_shift.sort(key=lambda elem: float(elem[1][sorting]))    
stations= [ts[0][0] for ts in tw_shift[::len(tw_range)]]

cm_colors = [(0.7,0,0), (1, 0, 0), (0.9,0.3,0.4), (1, 1, 1), (0.2,0.5,1), (0, 0, 1), (0.2,0,0.5)]  # R -> W -> B
cmap_name = 'rgb'
rwb = LinearSegmentedColormap.from_list(cmap_name, cm_colors, N=9)
ar=0.2

fig = plt.figure(figsize=(12,5))
for i in range(0,len(st_syn)): 
    ratio=[]
    for sta in stations:
        m = [ts for ts in tw_shift if ts[0][0]==sta]
        ratio.append([np.max(np.abs(ts[4]))/np.max(np.abs(ts[i+5])) for ts in m])

    ax = fig.add_subplot(3,4,i+1)
    im = ax.imshow(ratio, interpolation='none', aspect=ar, cmap=rwb, vmin=0, vmax=2)
    ax.set_xticks(np.arange(0, len(tw_range), 1)) # Major ticks
    ax.set_yticks(np.arange(0, len(stations), 1))
    ax.set_xticks(np.arange(-.5, len(tw_range), 1), minor=True) # Minor ticks
    ax.set_yticks(np.arange(-.5, len(stations), 1), minor=True)
    ax.set_xticklabels(tw_range, rotation='horizontal', fontsize=8) # Labels for major ticks
    ax.set_yticklabels(stations, rotation='horizontal', fontsize=7)
    ax.grid(which='minor', color='k', linestyle='-', linewidth=0.5) # Gridlines based on minor ticks
    ax.set_title('data/syn with %s'%syn_label[i])
    ax.set_xlabel('$tw_1$ (hrs)')

xlim=ax.get_xlim()
ax = fig.add_axes(ax.get_position(), frameon=False); ax.set_aspect(ar)
ax.tick_params(labelbottom='off',labeltop='off', labelleft="off", labelright='on',
               bottom='off', left='off', top='off', right='off', which='major')
ax.tick_params(labelbottom='off',labeltop='off', labelleft="off", labelright='on',
               bottom='off', left='off', top='off', right='on', which='minor')

ax.set_xlim(xlim)
ax.set_yticks(np.arange(0, len(stations), 1))
ax.set_yticks(np.arange(-.5, len(stations), 1), minor=True)
ax.set_yticklabels([int(ts[1][sorting]) for ts in tw_shift[0::len(tw_range)]], 
                    rotation='horizontal', fontsize=7)
ax.yaxis.set_label_position("right")
if sorting == 0:
    ax.set_ylabel('latitude', rotation=-90)
elif sorting == 1:
    ax.set_ylabel('longitude', rotation=-90)
elif sorting == 2:
    ax.set_ylabel('BAZ', rotation=-90)
elif sorting == 3:
    ax.set_ylabel('Epicentral distance', rotation=-90)
    
cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4]) 
clb = plt.colorbar(im, orientation='vertical', cax=cbar_ax)
tick_locator = ticker.MaxNLocator(nbins=10)
clb.locator = tick_locator
clb.update_ticks()
clb.ax.tick_params(labelsize=8) #set tick size
clb.ax.get_xaxis().set_ticks([]); y = [-0.1, 1.1]
for j, lab in enumerate(['overestimating Q','underestimating Q']):
    clb.ax.text(0.5, y[j], lab, ha='center', va='center')
fig.get_axes()[0].annotate('%s with $tw_2$=%s hrs, %s/%s stations'%(mode[0], tw_shift[0][2][1], len(stations), len (st)), 
                            (0.5, 0.9), xycoords='figure fraction', ha='center', fontsize=20) 
plt.show()
 
#-----  map plot of data/syn amplitude ratio
t=4
fig = plt.figure(figsize=(12,5))
for i in range(0,len(st_syn)): 
    ratio=[]
    for sta in stations:
        m = [ts for ts in tw_shift if ts[0][0]==sta]
        ratio.append([np.max(np.abs(ts[4]))/np.max(np.abs(ts[i+5])) for ts in m])

    ax = fig.add_subplot(3,4,i+1)
    m = Basemap(projection='moll',lon_0=180,resolution='c')    
    m.drawcoastlines()
    m.fillcontinents(color='lightgray',lake_color='whitesmoke')
    m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
    m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
    m.drawmapboundary(fill_color='whitesmoke')
    
    for j in range(len(stations)):
        latlon = [ts for ts in tw_shift if ts[0][0]==stations[j]]
        x, y = m(latlon[0][1][1], latlon[0][1][0])
        scatter = m.scatter(x, y, c=ratio[j][t], vmin=0, vmax=2,
                            s=200, marker="o", alpha=0.95, cmap=rwb, zorder=10)
    
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
    fig.add_axes(cax)
    clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
    clb.ax.set_title('data/syn ratio', position=(0.5, -2))
    ax.set_title('%s with tw:%s-%s'%(syn_label[i],tw_shift[t][2][0],tw_shift[t][2][1]), fontsize=18)
    plt.show()

#----- colormap misfit
sorting=0  #sorting=0  by latitude; sorting=2 by baz; sorting=3 by epicentral distance 
tw_shift.sort(key=lambda elem: float(elem[1][sorting]))    #
stations= [ts[0][0] for ts in tw_shift[::len(tw_range)]]

cm_colors = [(0.2,0,0.5), (0, 0, 1), (0.2,0.5,1), (1, 1, 1), (0.9,0.3,0.4), (1, 0, 0), (0.7,0,0)]  # B -> W -> R
cmap_name = 'bgr'
bwr = LinearSegmentedColormap.from_list(cmap_name, cm_colors, N=9)
ar=0.2
cmplx=0

if cmplx==0: mf_type='amplitude misfit'
if cmplx==1: mf_type='complex misfit'

fig = plt.figure(figsize=(12,5))
for i in range(0,len(st_syn)): 
    mf=[]
    for sta in stations:
        m = [1-ts[3][cmplx][i]/100 for ts in tw_shift if ts[0][0]==sta]
        mf.append([ts for ts in m])

    ax = fig.add_subplot(3,4,i+1)
    im = ax.imshow(mf, interpolation='none', aspect=ar, cmap=bwr, vmin=0, vmax=0.5)
    ax.set_xticks(np.arange(0, len(tw_range), 1)) # Major ticks
    ax.set_yticks(np.arange(0, len(stations), 1))
    ax.set_xticks(np.arange(-.5, len(tw_range), 1), minor=True) # Minor ticks
    ax.set_yticks(np.arange(-.5, len(stations), 1), minor=True)
    ax.set_xticklabels(tw_range, rotation='horizontal', fontsize=8) # Labels for major ticks
    ax.set_yticklabels(stations, rotation='horizontal', fontsize=7)
    ax.grid(which='minor', color='k', linestyle='-', linewidth=0.5) # Gridlines based on minor ticks
    ax.set_title('%s with %s'%(mf_type,syn_label[i]))
    ax.set_xlabel('$tw_1$ (hrs)')

xlim=ax.get_xlim()
ax = fig.add_axes(ax.get_position(), frameon=False); ax.set_aspect(ar)
ax.tick_params(labelbottom='off',labeltop='off', labelleft="off", labelright='on',
               bottom='off', left='off', top='off', right='off', which='major')
ax.tick_params(labelbottom='off',labeltop='off', labelleft="off", labelright='on',
               bottom='off', left='off', top='off', right='on', which='minor')

ax.set_xlim(xlim)
ax.set_yticks(np.arange(0, len(stations), 1))
ax.set_yticks(np.arange(-.5, len(stations), 1), minor=True)
ax.set_yticklabels([int(ts[1][sorting]) for ts in tw_shift[0::len(tw_range)]], 
                    rotation='horizontal', fontsize=7)
ax.yaxis.set_label_position("right")
if sorting == 0:
    ax.set_ylabel('latitude', rotation=-90)
elif sorting == 1:
    ax.set_ylabel('longitude', rotation=-90)
elif sorting == 2:
    ax.set_ylabel('BAZ', rotation=-90)
elif sorting == 3:
    ax.set_ylabel('Epicentral distance', rotation=-90)
    
cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4]) 
clb = plt.colorbar(im, orientation='vertical', cax=cbar_ax)
tick_locator = ticker.MaxNLocator(nbins=10)
clb.locator = tick_locator
clb.update_ticks()
clb.ax.tick_params(labelsize=8) #set tick size
clb.ax.get_xaxis().set_ticks([]); y = [-0.1, 1.1]
for j, lab in enumerate(['overestimating Q','underestimating Q']):
    clb.ax.text(0.5, y[j], lab, ha='center', va='center')
fig.get_axes()[0].annotate('%s with $tw_2$=%s hrs, %s/%s stations'%(mode[0], tw_shift[0][2][1], len(stations), len(st)), 
                            (0.5, 0.9), xycoords='figure fraction', ha='center', fontsize=20) 
plt.show()
 
#-----  map plot misfit
t=3 # selecting which tw2 to plot
cmplx=0
if cmplx==0: mf_type='amplitude misfit'
if cmplx==1: mf_type='complex misfit'

fig = plt.figure(figsize=(12,5))
for i in range(0,len(st_syn)): 
    mf=[]
    for sta in stations:
        m = [1-ts[3][cmplx][i]/100 for ts in tw_shift if ts[0][0]==sta]
        mf.append([ts for ts in m])

    ax = fig.add_subplot(3,4,i+1)
    m = Basemap(projection='moll',lon_0=180,resolution='c')    
    m.drawcoastlines()
    m.fillcontinents(color='lightgray',lake_color='whitesmoke')
    m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
    m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
    m.drawmapboundary(fill_color='whitesmoke')
    
    for j in range(len(stations)):
        latlon = [ts for ts in tw_shift if ts[0][0]==stations[j]]
        x, y = m(latlon[0][1][1], latlon[0][1][0])
        scatter = m.scatter(x, y, c=mf[j][t], vmin=0, vmax=1,
                            s=150, marker="o", alpha=0.95, cmap=bwr, zorder=10)
    
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
    fig.add_axes(cax)
    clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
    clb.ax.set_title('%s'%mf_type, position=(0.5, -2))
    ax.set_title('%s with tw:%s-%s'%(syn_label[i],tw_shift[t][2][0],tw_shift[t][2][1]), fontsize=18)
    plt.show()
