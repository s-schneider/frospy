#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:35:36 2019

@author: talavera
"""
import numpy as np
import obspy
from obspy import read as read
from frospy.util.array_util import stack
from frospy.core.segment import read as read_seg
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from obspy.core.utcdatetime import UTCDateTime

def func(x, A0, q):
#    f=0.814368*1e-3 # PREM
    f=0.814601*1e-3 # SC measurement
    om0 = f * 2 * np.pi * 3600 *24 # in rad/day
    return A0 * np.exp(-om0*x*0.5*q)

events = [
          "011307Z",  # aftershock 10 March, 2007, 2 stations
          "022710Z",  # good!
          "031111Z",  # fit witouth trimmning
          "032805Z",  # after removing stations good fit!
          "050306Z",  # remove? after shock, 4 July 2006?
#          "052413Z",  # not included in inversion
          "122604Z",  # fit w/o trimming, aftershock 15, march 2005?
          ]
home = "/net/home/talavera/eejit/splitting/00s00/Qcevents/self/c00=1"
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
for e in events:
    seg_file = "%s/%s/%s.dat"%(home,e,e)
    seg = read_seg(seg_file)
    stations = [s.stats.station for s in seg]
#    if e=="122604Z":
#         # Smaller ntps Remove from inversion!
#        stations.remove("STU")
#        stations.remove("VTS")
#        stations.remove("RGN")
#        # Okal & Stein, 2009 sued stations
#        # CTA, MAJO, PPT and YSS 
#        # From this I use: PPT, CTAO
#    if e=="032805Z":
#        # aftershock? Remove from inversion! done!
#        stations.remove("ANTO")
#        stations.remove("DWPF")
#        stations.remove("PAB")
#        stations.remove("ESK")
#    if e=="050306Z":
#        # aftershock? Remove from inversion! done!
#        stations.remove("ALE")
    st = read('/net/home/talavera/eejit/data/VHZ/%s.ahx'%e)
    st_filt = st.copy()
    
    tmax = st_filt[0].stats.npts*10/3600
#    tmax = 734400*10/3600
    start = st_filt[0].stats.starttime
    end = st_filt[0].stats.endtime
#    print(start,end)
    tw1 = seg[0].tw1
    tw2 = seg[0].tw2
    st_filt.trim(start+tw1*3600,end-(tmax-tw2)*3600)
   
    data_for_stack = [] 
    for tr in st_filt:
        station = tr.stats.station
        if station in stations:
            fw1 = seg.select(station=station)[0].fw1*1e-3
            fw2 = seg.select(station=station)[0].fw2*1e-3
#            tw1 = seg.select(station=station)[0].tw1
#            tw2 = seg.select(station=station)[0].tw2
#            tmax  = tr.stats.npts*10/3600
#            start = tr.stats.starttime
#            end   = tr.stats.endtime        
            tr.filter('bandpass', freqmin=fw1, freqmax=fw2, 
                      corners=2, zerophase=True)
            data_for_stack.append(tr.data)
#            tr.plot()
    
    data_for_stack = np.array(data_for_stack)
    
    # Stations have to have the same sampling rate!
    stacked_data = stack(data_for_stack)
    
    data_envelope = obspy.signal.filter.envelope(stacked_data)
    A0 = max(data_envelope)
    index = np.argmax(data_envelope)
    data_envelope=data_envelope[index::]
    
    # x is tranformed from seg to days and is reflected in func
    x=np.linspace(tw1/24,tw2/24,num=len(data_envelope))
#    ax1.plot(x,np.log10(data_envelope/max(data_envelope)),label=e)
    ax1.plot(x,data_envelope,label=e)
    ax2.plot(x,data_envelope)
    
    popt, pcov = curve_fit(func, x, data_envelope, p0=(A0, 1./5327.))
    ax2.plot(x, func(x, *popt), color="k")
    print("%s, Q=%2.f"%(e, 1/popt[1]))
#popt, pcov = curve_fit(func, x, np.log10(data_envelope/max(data_envelope)))
ax1.plot(x, func(x, *popt), color="k", label="regression")
ax1.set_yscale('log')
#ax.set_xlim(-1,85)
ax1.legend(ncol=1,loc="lower left")
plt.show()

#011307Z, Q=6802
#022710Z, Q=5976
#031111Z, Q=6063
#032805Z, Q=6090
#122604Z, Q=6005
#050306Z, Q=7916

events = [
          "011307Z",  # aftershock 10 March, 2007, 2 stations
          "022710Z",  # good!
          "031111Z",  # fit witouth trimmning
          "032805Z",  # after removing stations good fit!
          "050306Z",  # remove? after shock, 4 July 2006?
          "052413Z",  # not included in inversion
          "122604Z",  # fit w/o trimming, aftershock 15, march 2005?
          ]
for e in events:
    print(e)
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    seg_file = "%s/%s/%s.dat"%(home,e,e)
    seg = read_seg(seg_file)
    stations = [s.stats.station for s in seg]
    for sta in stations:
        st = read('/net/home/talavera/eejit/data/VHZ/%s.ahx'%e)
        st_filt = st.copy()
        
        tmax = st_filt[0].stats.npts*10/3600
    #    tmax = 734400*10/3600
        start = st_filt[0].stats.starttime
        end = st_filt[0].stats.endtime
    #    print(start,end)
        tw1 = seg[0].tw1
        tw2 = seg[0].tw2
        st_filt.trim(start+tw1*3600,end-(tmax-tw2)*3600)
       
        fw1 = seg.select(station=sta)[0].fw1*1e-3
        fw2 = seg.select(station=sta)[0].fw2*1e-3   
        tr = st_filt.select(station=sta)
        # Maybe is better to first filter and trim after filtering
        tr.filter('bandpass', freqmin=fw1, freqmax=fw2, 
                  corners=2, zerophase=True)
        data=tr[0].data
    #   tr[0].plot()

        
        data_envelope = obspy.signal.filter.envelope(data)
        A0 = max(data_envelope)
        index = np.argmax(data_envelope)
        data_envelope=data_envelope[index::]
        
        # x is tranformed from seg to days and is reflected in func
        x=np.linspace(tw1/24,tw2/24,num=len(data_envelope))
    #    ax1.plot(x,np.log10(data_envelope/max(data_envelope)),label=e)
        
        popt, pcov = curve_fit(func, x, data_envelope, p0=(A0, 1./5327.))
        ax2.plot(x, func(x, *popt), color="k")
        print("%s %2.f"%(sta, 1/popt[1]))
        
        ax1.plot(x,data_envelope,label="%s: Q=%2.f"%(sta,1/popt[1]))
        ax2.plot(x,data_envelope)
    #popt, pcov = curve_fit(func, x, np.log10(data_envelope/max(data_envelope)))
    ax1.plot(x, func(x, *popt), color="k", label="regression")
    ax1.set_yscale('log')
    #ax.set_xlim(-1,85)
    ax1.legend(ncol=1,loc="lower left")
    plt.suptitle(e)
    plt.show()