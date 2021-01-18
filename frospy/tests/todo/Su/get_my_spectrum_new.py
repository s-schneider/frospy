#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:05:55 2019

@author: talavera
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:26:29 2018

@author: talavera
"""
import numpy as np
from obspy.core import AttribDict
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from frospy.util.read import read as read_modes
from frospy.spectrum.spectrum_interactive import set_defaults_init_main
from frospy.spectrum.app import check_input
from frospy.spectrum.spectrum_interactive import set_defaults_init_main
from frospy.spectrum.controllers import (read_data, read_syn, load_cmt,
                                       _plot_modes_in_ax, get_mode,
                                       set_plot_args, process_input, 
                                       search_stream, report_input)
from frospy.util.base import (fourier_transform, get_times, taper)
from frospy.core.segment import read as read_segments
from frospy.util.base import inv4stream
from obspy import Stream
from frospy.core import Spectrum
from frospy.core.spectrum.plot import plot_spectrum

event='031111B'
data='PREM/%s.ahx'%event
syn=[
     'PREM/%s.ahx.syn'%event, # self
     'HT/%s.ahx.syn'%event,
     ] 

syn_label=[
           'PREM',
           'He & Tromp, 1996', 
           'My Inversion'
           ]

segment = read_segments('segments/%s.dat'%event)
stations = [s.station for s in segment]

#mode=['${}_1 S_0$']; fw = [1.618, 1.641]; tw = [5,280]
#mode=['${}_2 S_0$']; fw = [2.489, 2.534]; tw = [30,150]
#mode=['${}_3 S_0$']; fw = [3.255, 3.285]; tw = [35,105]
mode=['${}_4 S_0$']; fw = [4.085, 4.125]; tw = [25,100]
#mode=['${}_5 S_0$']; fw = [4.872, 4.895]; tw = [35,90]
#mode=['${}_6 S_0$']; fw = [5.724, 5.759]; tw = [20,75]


iargs = AttribDict()
iargs.syn = syn
iargs.syn_label = syn_label
iargs.fw = [4.085, 4.125]
iargs.tw = [25,100]
iargs.station = "ANMO"

data, iargs = check_input(data, iargs)

main = AttribDict(iargs)
main = set_defaults_init_main(main)

if main.cmt_file:
    main.cmt = load_cmt(file=main.cmt_file, verbose=main.verbose)
else:
    main.cmt = load_cmt(id=main.cmt_id, verbose=main.verbose)

# Loading pick file
if main.segment_file is not None:
    main.segments = load_segments([main.segment_file])
else:
    main.segments = None
    
main = read_data(data, main)
main.inv = inv4stream(main.st)

report_input(main)

main = search_stream(main, main.station)
main.tr = main.st_work[0]
# set window parameters ################################
if main.segments is not None:
    try:
        stat = main.tr.stats.station
        pick = main.segments.select(station=stat)[main.j]
        if main.set_seg_tw:
            main.tw = [pick.tw1, pick.tw2]
        if main.set_seg_fw:
            main.fw = [pick.fw1 - 0.05, pick.fw2 + 0.05]
    except IndexError:
        msg = 'No picked segments in file for %s\n'
        print(msg % main.tr.stats.station)
        pass

main.syn_trs = Stream()
if main.syn:
    for s in main.st_syn_work:
        main.syn_trs += s.select(station=main.tr.stats.station)

spec = Spectrum(main.tr, main.tw[0], main.tw[1], main.syn_trs,
                taper=main.taper_shape,
                syn_label=main.syn_label, label=main.data_label)

main.startlabel = spec.flabel(main.fw[0])
main.endlabel = spec.flabel(main.fw[1])

plot_args = set_plot_args(main, spec)
plotargs = plot_spectrum(plot_args, main.i)