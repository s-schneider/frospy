#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:26:29 2018

@author: talavera
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from frospy.util.base import format_exponent
from frospy.util.read import read_modes
from frospy.preprocessing.spectrum import read_data, read_syn, get_times,taper, fourier_transform
from frospy.preprocessing.spectrum import _plot_modes_in_ax, get_mode

def triang(t):
    if t >= -1.0 and t < 0.0:
        return 1.0 + t
    elif t >= 0.0 and t <= 1.0:
        return 1.0 - t
    return 0.0

def integral_triang(t):
    if t <= -1.0:
        return 0
    elif t >= -1.0 and t < 0.0:
        return 0.5 * (t + 1)**2
    elif t >= 0.0 and t <= 1.0:
        return -0.5 * t**2 + t + 0.5
    else:
        return 1.0

def my_fourier_transform(trace, Htstart, Htend, nstart, nend, shape, mfac=1):
    maxdata = 65536  # len(trace.data)  # 65536  why is it this value?
    indvec = np.arange(maxdata)
    delomeg = 2. * np.pi / (maxdata * trace.stats.delta)
    f = 1000. * indvec * delomeg / (2. * np.pi)
    f = f - (1000. * delomeg / (2. * np.pi))
    mid = Htstart+0.5*(Htend-Htstart)
    cfac = np.exp(1j * f * mid * 2. * np.pi/1000.)
    # no filter
    #tracefill = taper(nstart, nend, trace, shape) * mfac
    X = np.fft.fft(trace, maxdata)
    FT = X*cfac
    phase = np.angle(FT)

    return f, FT, phase, delomeg


event='031111B'
data='/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx'%event
syn=['/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx.syn'%event] #s20rts
#'/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/prem/%s.ahx.syn'%event]
syn_label=['0-7.5mHz'] #-S20RTS', '0-7.5mHz-PREM']
fw=[3.9, 4.2]; tw=[10,50] #f_4s0=4106.20
#halfDuration=70.0; 
halfDuration=(1.05e-08 * (5.31e+29)**(1./3))
duration=halfDuration*2.0
amp = 1.0#/halfDuration

taper_shape='hanning'
Htstart = tw[0] * 3600.
Htend = tw[1] * 3600.
Htmid = (tw[1] - tw[0])/2.
wstart = fw[0]
wend = fw[1]

Htstart_org, Htend_org, Htmid_org = Htstart, Htend, Htmid

st, st_work, syn = read_data(data, syn)
st_syn, st_syn_work = read_syn(syn, st, syn_label)

#for trace in tr:
#for j in range(len(tr)):
#for i in range(2):
i = 0; key = ['nps']
while key != 'quit':
    # setting-up data
    tr_data = st_work[i]
    tr_syn = st_syn_work[0][i]
    times = get_times(tr_data, Htstart, Htend, Htstart_org, Htend_org)
    tstart, nstart, nend, t, Htstart, Htend, Htmid, Terr = times

    # filtering sesimograms in time
    trfill = taper(nstart, nend, tr_data, taper_shape)
    trsynfill = taper(nstart, nend, tr_syn, taper_shape)
    
    Fxx = []; Fxxsyn = []; Fxxsyn_stf = []
    xp = []; xpsyn = []; xpsyn_stf = []
    startlabel = []
    endlabel = []
    freq = []

    # Source Time Function (STF)
    dt=tr_data.stats.delta
    #nt1 = int(np.ceil(duration) / dt) + 1 
    nt1 = tr_data.stats.npts
    #stf = [amp * triang((j*dt) / halfDuration - 1.0) for j in range(0, nt1)]
    stf = [amp * integral_triang((j*dt) / halfDuration - 1.0) for j in range(0, nt1)]
    #convol = np.convolve(tr_syn, stf, mode='same') # CONVOLVING STF with Syn
    #syn_stf = tr_data.copy(); syn_stf.data = convol
    #convolfill = taper(nstart, nend, syn_stf, taper_shape)
    syn_stf = tr_data.copy(); syn_stf.data = np.array(stf)
    
    # Perform Fouriertransform and filter data
    f, FT, XP, delomeg = fourier_transform(tr_data, Htstart, Htend, nstart, nend, taper_shape, 1)
    freq.append(f)
    Fxx.append(FT)
    xp.append(XP)
    
    # Perform Fouriertransform syn
    void, FT, XP, void = fourier_transform(tr_syn, Htstart, Htend, nstart, nend, taper_shape, 1)
    Fxxsyn.append(FT)
    xpsyn.append(XP)
    
    # Perform Fouriertransform convolution STF with Syn
    void, FT, XP, void = fourier_transform(syn_stf, Htstart, Htend, nstart, nend, taper_shape, 1)
    Fxxsyn_stf.append(FT)
    xpsyn_stf.append(XP)
    #Fxxsyn_stf = [Fxxsyn_stf[0]*Fxxsyn[0]]
    
    # setting limits
    startlabel.append(int(wstart * 2. * np.pi/(1000. * delomeg)) + 1)
    endlabel.append(int(wend * 2. * np.pi/(1000. * delomeg)) + 1)
    startlabel=startlabel[0]
    endlabel=endlabel[0]
    
    #plotting
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(3,1,1)
    ax.plot(f[startlabel:endlabel+1], xp[0][startlabel:endlabel+1],
            linestyle='solid', color='black', linewidth=.75)
    ax.plot(f[startlabel:endlabel+1], xpsyn[0][startlabel:endlabel+1],
            linestyle='dashed', color='blue', linewidth=.75)
    ax.plot(f[startlabel:endlabel+1], xpsyn_stf[0][startlabel:endlabel+1],
            linestyle='dashed', color='red', linewidth=.75)
    box = ax.get_position() # Shrink current axis's height by 35% on the bottom
    ax.set_position([box.x0, box.y0 + box.height * 0.65,
                         box.width, box.height * 0.35])
    ax.set_ylabel('Phase')
    ax.get_xaxis().set_ticklabels([])
    ax.set_xlim(f[startlabel], f[endlabel+1])
    ax.set_title(r'%s, tw: %s - %s' % (event, tw[0], tw[1]), fontsize=15)

    ax = fig.add_subplot(3,1,2)
    dlabel = '%s' % (tr_data.stats.station)
    ax.plot(f[startlabel:endlabel+1], abs(Fxx[0][startlabel:endlabel+1]),
            linestyle='solid', color='black', linewidth=0.75, label=dlabel)
    ax.plot(f[startlabel:endlabel+1], abs(Fxxsyn[0][startlabel:endlabel+1]),
            linestyle='dashed', color='blue', linewidth=0.75, label=dlabel+'-syn')
    ax.plot(f[startlabel:endlabel+1], abs(Fxxsyn_stf[0][startlabel:endlabel+1]),
            linestyle='dashed', color='red', linewidth=0.75, label=dlabel+'-conv')
    box = ax.get_position() # Shrink current axis's height by 10% on the bottom
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 1.5])
    modedata = read_modes()#; modes = modedata
    modes=get_mode(['4s0', '11S2', '10S2'], modedata)
    f_lim = [f[startlabel], f[endlabel+1]]
    _plot_modes_in_ax(ax, modes, f_lim)
    ax.set_ylabel('Amplitude')
    ax.set_xlabel('frequency (mHz)')
    ax.legend()
    mf = ticker.ScalarFormatter(useMathText=True)
    mf.set_powerlimits((-2,2)) # for exponents greater than +-2
    plt.gca().yaxis.set_major_formatter(mf)
    #ax = format_exponent(ax)

    ax = fig.add_subplot(3,1,3)
    ax.plot(t[nstart:nend] / 3600., trfill[nstart:nend], linestyle='solid', color='black', linewidth=0.75, label=dlabel)
    #ax.plot(t[nstart:nend] / 3600., convolfill[nstart:nend], linestyle='dashed', color='red', linewidth=0.75, label=dlabel+'-conv')
    ax.plot(t[nstart:nend] / 3600., trsynfill[nstart:nend], linestyle='dashed', color='blue', linewidth=0.75, label=dlabel+'-syn')
    box = ax.get_position() # Shrink current axis's height by 10% on the bottom
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
                     box.width, box.height * 0.7])
    ax.set_ylabel('Amplitude')
    ax.set_xlabel('time (hrs)')
    mf = ticker.ScalarFormatter(useMathText=True)
    mf.set_powerlimits((-2,2))
    plt.gca().yaxis.set_major_formatter(mf)
    plt.show()

    key = raw_input('Please an action (nsp, bts, psp, tw, fw, quit): ')
    if i >= len(st_work): i = 0
    elif key == 'next' or key == 'nsp' or key == "": i += 1
    elif key == 'quit' or key == 'exit' or key == 'q': break
    elif key == 'bts': i = 0
    elif key == 'fw':
        fw1 = raw_input('Input start fw: '); wstart=float(fw1)
        fw2= raw_input('Input end fw: '); wend=float(fw2)
    elif key == 'tw':
        tw1 = raw_input('Input start tw: '); tw[0]=float(tw1)
        tw2 = raw_input('Input end tw: '); tw[1]=float(tw2)
        Htstart = tw[0] * 3600.
        Htend = tw[1] * 3600.
        Htmid = (tw[1] - tw[0])/2.
    elif key == 'psp':
        if i != 0:
            i -= 1
        else:
            i = 0
            print('Already at the beginning.\n')
    #ax = format_exponent(ax)
    #plt.savefig('test', transparent=True, bbox_inches='tight', pad_inches=0, dpi=300)
    

#from frospy.util.read import read_st, read_modes
#import matplotlib.ticker as ticker
#from obspy.core import AttribDict
#    
#def read_data(data, syn):
#    st = []
#    if type(data) is list:
#        syn = None
#        for file in data:
#            st_tmp = read_st(file, 'ah')
#            st.append(st_tmp)
#        st_work = st[:]
#
#    else:
#        st_tmp = read_st(data, 'ah')
#        st_tmp.sort(['station'])
#        st = st_tmp
#        st_work = st.copy()
#    print("\t data loaded")
#
#    return st, st_work, syn
#
#def read_syn(syn, st, syn_label):
#    st_syn = []
#    if type(syn) is str:
#        stsyn_tmp = read_st(syn)
#        stsyn_tmp.sort(['station'])
#        st_syn.append(stsyn_tmp)
#        print("\t (%i/%i) synthetics loaded" % (1, len(st_syn)))
#
#    elif type(syn) is list:
#        for _i, file in enumerate(syn):
#            stsyn_tmp = read_st(file)
#            st_syn.append(stsyn_tmp)
#            print("\t (%i/%i) synthetics loaded" % (_i+1, len(syn)))
#
#    elif syn is None:
#        return None, None
#
#    else:
#        msg = "Wrong input format of synthetic data files."
#        msg += "Must be string or list of strings"
#        raise IOError(msg)
#
#    for syns in st_syn:
#        if len(syns) != len(st):
#            msg = 'Number of synthetic traces differs from data traces \n'
#            raise IOError(msg)
#        elif syn_label is not None and len(syn) != len(syn_label):
#            msg = 'Number of syn_labels (%i) differs from number of '
#            msg += 'synthetic traces (%i) \n'
#            msg = msg % (len(syn_label), len(syn))
#            raise IOError(msg)
#    st_syn_work = st_syn[:]
#
#    return st_syn, st_syn_work
#
#def starttime(trace):
#    rsec = trace.stats.starttime.second
#    rsec = rsec + trace.stats.starttime.microsecond * 1E-6
#    esec = trace.stats.ah.event.origin_time.second
#    esec = esec + trace.stats.ah.event.origin_time.microsecond * 1E-6
#    rmn = trace.stats.starttime.minute
#    emn = trace.stats.ah.event.origin_time.minute
#    rhr = trace.stats.starttime.hour
#    ehr = trace.stats.ah.event.origin_time.hour
#    nde = trace.stats.starttime.day
#    ndr = trace.stats.ah.event.origin_time.day
#
#    tstart = rsec-esec+60.*((rmn-emn) + 60.*((rhr-ehr) + 24.*(nde-ndr)))
#
#    return tstart
#
#def get_times(tr, Htstart, Htend, Htstart_org, Htend_org):
#
#    # Set errors to false
#    STshort = False
#    STlong = False
#    ETlong = False
#    err = False
#    tstart = starttime(tr)
#
#    if Htstart_org is not None:
#        Htstart, Htend = Htstart_org, Htend_org
#
#    Htmid = (Htstart + Htend) / 2.
#    nstart = int(round((Htstart-tstart)/tr.stats.delta)+1)
#    nend = int(round((Htend-tstart)/tr.stats.delta))
#    t = np.arange(tr.stats.npts) * tr.stats.delta
#    t = t + tstart
#
#    if nstart < 0:
#        Tdiff = Htend - Htstart
#        Htstart = tstart
#        Htend = Htstart + Tdiff
#        Htmid = Htstart + Tdiff/2.
#        nstart = 0
#        nend = int(round((Htend-tstart)/tr.stats.delta))
#        t = t + tstart
#        s = Htstart / 3600.
#        e = Htend / 3600.
#        STshort = True
#
#    elif nstart > len(tr.data):
#        Tdiff = tr.stats.endtime - tr.stats.starttime
#        Htstart = tr.stats.starttime - tr.stats.ah.event.origin_time
#        Htend = Htstart + Tdiff
#        Htmid = Htstart + Tdiff/2.
#        nstart = 0
#        nend = int(round((Htend-tstart)/tr.stats.delta))
#        s = Htstart / 3600.
#        e = Htend / 3600.
#        STlong = True
#
#    if nend > tr.stats.npts:
#        nend = tr.stats.npts - 1
#        Htend = np.floor(nend * tr.stats.delta + tstart)
#        Tdiff = Htend - Htstart
#        Htmid = Htstart + Tdiff/2.
#        ETlong = True
#
#    s = Htstart / 3600.
#    e = Htend / 3600.
#
#    msg = ''
#    if STshort or STlong or ETlong:
#        err = True
#
#    if STshort:
#        msg += '\n\033[93mStarttime earlier then record time'
#
#    if STlong:
#        msg += '\n\033[93mStarttime later then record time'
#
#    if ETlong:
#        msg += '\n\033[93mEndtime later then record time'
#    if msg:
#        msg += '\nTimewindow set to %i-%i h\033[0m\n'
#        print(msg % (s, e))
#
#    return(tstart, nstart, nend, t, Htstart, Htend, Htmid, err)
#
#def taper(nstart, nend, trace, shape='hanning'):
#    npun = nend-nstart
#    # ppi = 2*np.pi/(npun-1)
#    # ni = npun/2
#    maxdata = len(trace.data)
#    if shape.lower() == 'hanning':
#        # fillpart = 0.5 + 0.5 * np.cos(ppi * (np.arange(1, npun)-ni))
#        fillpart = np.hanning(npun)
#    elif shape.lower() == 'bartlett':
#        fillpart = np.bartlett(npun)
#    elif shape.lower() == 'hamming':
#        fillpart = np.hamming(npun)
#    elif shape.lower() == 'blackman':
#        fillpart = np.blackman(npun)
#    elif shape.lower() == 'kaiser':
#        fillpart = np.kaiser(npun, 5)
#    elif shape.lower() == 'boxcar':
#        fillpart = np.ones(npun)
#
#    fill = np.zeros(maxdata)
#    fill[nstart:nend] = fillpart
#    tracefill = trace.data * fill
#
#    return tracefill
#
#def fourier_transform(trace, Htstart, Htend, nstart, nend, shape, mfac=1):
#    maxdata = 65536  # len(trace.data)  # 65536  why is it this value?
#    indvec = np.arange(maxdata)
#    delomeg = 2. * np.pi / (maxdata * trace.stats.delta)
#    f = 1000. * indvec * delomeg / (2. * np.pi)
#    f = f - (1000. * delomeg / (2. * np.pi))
#    mid = Htstart+0.5*(Htend-Htstart)
#    cfac = np.exp(1j * f * mid * 2. * np.pi/1000.)
#    tracefill = taper(nstart, nend, trace, shape) * mfac
#    X = np.fft.fft(tracefill, maxdata)
#    FT = X*cfac
#    phase = np.angle(FT)
#
#    return f, FT, phase, delomeg
#
#def _plot_modes_in_ax(ax, modes, f_axis):
#    modes_sorted = sorted(modes.items(), key=lambda x: x[1]['freq'])
#    trans = ax.get_xaxis_transform()
#    for m in modes_sorted:
#        mode = m[1]
#        if mode.freq > f_axis[0] and mode.freq < f_axis[1]:
#            ax.text(mode.freq, 1.03, r"$%s$" % mode.name,
#                    transform=trans, rotation=45, ha="center",
#                    va='center')
#            ax.axvline(x=mode.freq, linestyle=':', linewidth=1, color='grey')
#    return
#
#def get_mode(names, modedata):
#    if names is None:
#        names = raw_input('Which mode? (name, e.g. 0S15) -->  ')
#    if type(names) is str:
#        names = names.split()
#
#    modes = AttribDict()
#    for name in names:
#        n = list(name[0:5])
#        # Checking for right format: Digit Letter Digit
#        if n[0] == '0' and n[1].isdigit():
#            n = n[1:]
#            n[1] = n[1].capitalize()
#
#        elif n[0].isdigit() and not n[1].isdigit():
#            n[1] = n[1].capitalize()
#
#        mode = str().join(n)
#        modes[mode] = modedata[mode].copy()
#
#    return modes