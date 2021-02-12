#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

from frospy.spectrum.controllers import (
                    load_cmt, read_data, load_segments, report_input,
                    get_spectrum, auto_picker, get_pick, get_misfit_snr,
                    set_plot_args
                    )

try:
    from frospy.core.spectrum.plot import plot_spectrum, plot_segments
except ImportError:
    pass
from frospy.core.modes import read as read_modes
from frospy.core.segment import Segment
from frospy.util.base import inv4stream
from frospy.plot.nmplt import mode_lines

import matplotlib.pyplot as plt

from obspy.core import AttribDict

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass


def run(data=None, **specargs):
    main = AttribDict(specargs)
    main = set_defaults_init_main(main)

    if data is not None:
        """
        If called using `frospy.preprocessing.spectrum.spectrum` with `runall`
        """
        # Loading cmt file
        if main.cmt_file:
            main.cmt = load_cmt(file=main.cmt_file, verbose=main.verbose)
        else:
            main.cmt = load_cmt(id=main.cmt_id, verbose=main.verbose)

        # Loading pick file
        if main.segment_file is not None:
            main.segments = load_segments([main.segment_file], main.verbose)

        pick = None

        # Don't check for identical data files but check later for correct
        # channel selection
        main = read_data(data, main)
        main.inv = inv4stream(main.st)

        if main.verbose:
            print(report_input(main))

    else:
        """
        If called internally from
        `frospy.preprocessing.spectrum._spectrum_interactive`
        all arguments are parsed as the dict `specargs`
        """
        main = specargs['main']
        pick = specargs['pick']

    # Convert from hours to seconds
    all_modes = read_modes()

    if main.verbose is True:
        delmsg = "%s\t| %s \033[91mdeleted\033[0m"
        savemsg = "%s\t| %s \033[92mkept\033[0m"
        print('Station | Misfit:\n ')

    if main.plot is True:
        specs = []
    for tr_i, tr in enumerate(main.st_work):
        """
        Picking behavior for autopick is set here.
        If modes is not None, fw will be ignored and the freq window is set
        according to the auto_picker.
        If modes is None the freq window is set according to fw[0] and fw[1]
        If segments is given as input the fw and tw values are defined by it
        """
        if specargs['segment_file'] is not None:
            pick = main.segments.select(station=tr.stats.station)[0]
            if main.set_seg_fw:
                fw = [pick.fw1 - 0.05, pick.fw2 + 0.05]
            else:
                fw = main.fw
            if main.set_seg_tw:
                tw_hour = [pick.tw1, pick.tw2]
            else:
                tw_hour = main.tw
        else:
            tw_hour = main.tw
            fw = main.fw

        if 'update_segments' in specargs:
            fw = main.fw

        main.spec = get_spectrum(tr, main.st_syn_work, tw_hour, tr_i,
                                 specargs['taper_shape'],
                                 specargs['stream_order'])

        if main.plot is True:
            specs.append(main.spec)
            continue

        if specargs['autopick'] is not None:
            fw_ap = auto_picker(main.spec, specargs['modes'])
        else:
            fw_ap = fw
        if fw_ap is None:
            continue
        pick = get_pick(main.spec, fw_ap, all_modes, specargs['weighting'],
                        event=main.cmt_id)

        # Calculate misfit
        fwhm, cmf, rmf, = get_misfit_snr(main.spec, pick)
        if cmf is not None:
            pick.stats.misfit = cmf[0]
        max_mf = specargs['max_mf']
        min_snr = specargs['min_snr']
        if max_mf is not None or min_snr is not None:

            if min_snr is not None:
                if pick.stats.snr < min_snr:
                    if specargs['segment_file'] is not None:
                        seg = main.segments
                        pick = seg.select(station=tr.stats.station)[0]
                        main.segments.remove(pick)
                        if main.verbose is True:
                            msg = delmsg % (main.spec.stats.station.code,
                                            cmf[0])
                            print(msg)
                    continue

            if max_mf is not None and cmf is not None:
                if cmf[0] > max_mf:
                    if specargs['segment_file'] is not None:
                        seg = main.segments
                        pick = seg.select(station=tr.stats.station)[0]
                        main.segments.remove(pick)
                        if main.verbose is True:
                            msg = delmsg % (main.spec.stats.station.code,
                                            cmf[0])
                            print(msg)
                    continue
            if main.verbose is True:
                msg = savemsg % (main.spec.stats.station.code, cmf[0])
                print(msg)

        # If there was no input, the segment file is created from scratch
        if main.segments is None:
            main.segments = Segment()

        if specargs['segment_file'] is None:
            main.segments += pick
        else:
            pick_old = main.segments.select(station=tr.stats.station)[0]
            main.segments.remove(pick_old)
            main.segments += pick

        if main.plot is True:
            main = plot_spectrum(main)

    if main.plot:
        fig, ax = plt.subplots()
        fig.set_size_inches(main.fig_size)
        ax.yaxis.set_ticks([0])
        for spec in specs:
            spec.plot(fw[0], fw[1], part='Amplitude', xlabel='frequency (mHz)',
                      ax=ax, ylabel='Amplitude', normalize=main.normalize)

        f = spec.stats.freq
        startlabel = spec.flabel(fw[0])
        endlabel = spec.flabel(fw[1])

        if main.modes is not None:
            mode_lines(ax, f[startlabel:endlabel+1], all_modes,
                       label_height=0.025,
                       label_width=0.03)
        if main.segments is not None:
            plot_segments(specs[0], main.segments, fw[0], fw[1], [ax])

        return fig, ax, main.segments

    if main.verbose:
        if main.segments is not None:
            print('No. of Segments: %s' % len(main.segments))

    return main.segments


def set_defaults_init_main(main):
    main.verbose = False
    main.recalc = True
    main.weighting = 'integrate'
    main.update_maps = True
    main.zoom_sl = 5.
    main.save_stream = False
    main.fs = 10.
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
    main.fs = 10

    try:
        main.normalize
    except Exception:
        main.normalize = False
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

    if main.seg_reference is True:
        main.set_seg_fw = True
        main.set_seg_tw = True
        if main.fw:
            main.set_seg_fw = False
    else:
        main.set_seg_fw = False
        main.set_seg_tw = True
        if main.fw:
            main.set_seg_tw = False

    main.i = 0
    main.i_old = None

    main.j = 0  # segment iterator
    if 'plot' not in main:
        main.plot = False

    if not main.fig_size:
        main.fig_size = (14, 7)

    return main
