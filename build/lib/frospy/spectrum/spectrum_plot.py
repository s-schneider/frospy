#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

from .controllers import (
                    load_cmt, read_data, load_segments, search_stream,
                    report_input, set_plot_args
                    )

from frospy.core import Spectrum
from frospy.core.modes import read as read_modes
from frospy.util.base import inv4stream
try:
    from frospy.core.spectrum.plot import plot_spectrum, plot_spectrum_partials
    NOWINDOW = False
except ImportError:
    NOWINDOW = True

from obspy.core import AttribDict
from obspy import Stream

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass


def run(data, **specargs):

    # create all local variables
    main = AttribDict(specargs)
    main = set_defaults_init_main(main)
    for key, value in specargs.items():
        main.update({key: value})

    print('Reading files ...')
    # Loading cmt file
    if main.cmt_file:
        main.cmt = load_cmt(file=main.cmt_file, verbose=main.verbose)
    else:
        main.cmt = load_cmt(id=main.cmt_id, verbose=main.verbose)

    # Loading pick file
    if main.segment_file is not None:
        main.segments = load_segments([main.segment_file])
    else:
        main.segments = None
    pick = None

    # Don't check for identical data files but check later for correct channel
    # selection
    main = read_data(data, main)
    main.inv = inv4stream(main.st)

    print(report_input(main))

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

    main.spec = Spectrum(main.tr, main.tw[0], main.tw[1], main.syn_trs,
                         taper=main.taper_shape,
                         syn_label=main.syn_label, label=main.data_label)

    main.startlabel = main.spec.flabel(main.fw[0])
    main.endlabel = main.spec.flabel(main.fw[1])

    # plot_args = set_plot_args(main, spec)
    if type(main.minispec) == list:
        main = plot_spectrum_partials(main)
    else:
        main = plot_spectrum(main)

    return main.rfig, main.rax  # plotargs.fig, plotargs.ax_old


def set_defaults_init_main(main):
    main.recalc = True
    main.weighting = 'integrate'
    main.update_maps = True
    main.zoom_sl = 5.
    main.save_stream = False
    main.noisewin = False
    main.overlap = False

    main.rmf = None
    main.cmf = None

    # Magnification of synthetics
    main.mfac = 1

    # Set reference figure and ax to None
    main.rfig = None
    main.rax = None
    main.seg_ax = None
    main.savefig = False
    try:
        main.fs
    except AttributeError:
        main.fs = 10
    main.border_width = 1
    main.tick_width = 1

    # Reading modes
    main.modes_all = read_modes()

    # Set Referece times
    if main.tw:
        main.tw_org = main.tw[:]
    else:
        main.tw_org = None
    main.qwindow = None
    main.min_snr = None
    main.max_mf = None

    if main.seg_reference:
        main.set_seg_fw = True
        main.set_seg_tw = True
        if main.fw:
            main.set_seg_fw = False
    else:
        main.set_seg_fw = False
        main.set_seg_tw = True

    main.i = 0
    main.i_old = None

    main.j = 0  # segment iterator

    return main
