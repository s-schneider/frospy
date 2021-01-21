#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os

from frospy.spectrum import spectrum_loop
from frospy.spectrum import spectrum_plot
from frospy.spectrum import spectrum_commandline
from frospy.spectrum import spectrum_gui

import obspy
try:
    from matplotlib.pyplot import cm
    NOWINDOW = False
except ImportError:
    NOWINDOW = True

from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass

"""
Module containing the functionality of the matlab scripts. Description
and documentation are about to follow.

Color Codes:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

Simon Schneider, 2017
"""


def spectrum(data=None, **iargs):
    """
    Program to work on spectra of seismic data

    Options for minispec:
     - True
     - 'amp'
     - 'lat'
     - list containing one or more of:
        'amp' : plot the amplitude spectrum
        'phase' : plot the phase
        'real' : plot real part of amplitude spectrum
        'imag' : plot Imaginary part of amplitude spectrum
        'gcp' : plot map with great circle path of station and event
        'stat' : plot map with all stations of file
        'lat' : plot amplitude spectrum against latitude
        'seis' : plot time series of data
        'mag' : plot cmt information of event


    Input:
    data=None,
    syn=None,
    syn_label=None,
    modes=False,
    data_label=None,
    segment_file=None,
    seg_reference=False,
    cmt_file=None,
    taper_shape='hanning',
    weighting='integrate',
    nowindow=False,
    minispec=False,
    fw=None,
    tw=None,
    min_snr=None,
    max_mf=None,
    search=None,
    cmap='rainbow',
    verbose=False,
    runall=False,
    autopick=None
    output=None

    Output:
    segments

    Example:
    data = '/data/simons/3-comp-data/120hrs/VHT/060994A.ahx'
    syn = []
    label = []
    for i in range(1,10):
        syn.append('Run_4/Iter_%i/T/synseis/060994A.ahx.syn' % i)
        label.append('Iter_%i' % i)

    segment_file = 'Run_4/Iter_1/T/synseis/060994A.dat'
    spectrum(data, syn, fw=[2.25, 2.42], tw=[3, 37], segment_file=segment_file,
             syn_label=label)
    """

    data, iargs = check_input(data, iargs)

    if type(iargs['search']) == str:
        iargs['search'] = [iargs['search']]

    if type(iargs['station']) is str:
        fig, ax = spectrum_plot.run(data, **iargs)
        return fig, ax

    elif iargs['runall'] is not False:
        out = spectrum_loop.run(data, **iargs)
        if 'plot' in iargs:
            fig, ax, segments = out[:]
            return fig, ax, segments
        else:
            segments = out
            return segments

    elif 'gui' in iargs:
        del iargs['gui']
        spectrum_gui.run(data, **iargs)

    else:
        segments = spectrum_commandline.run(data, **iargs)
        return segments


def get_defaults(args):

    defaults = {'syn': None,
                'syn_label': None,
                'modes': False,
                'segment_file': None,
                'seg_reference': False,
                'cmt_file': None,
                'taper_shape': 'hanning',
                'weighting': 'integrate',
                'nowindow': False,
                'minispec': False,
                'fw': None,
                'tw': None,
                'min_snr': None,
                'max_mf': None,
                'search': None,
                'cmap': 'rainbow',
                'cmap_highlight': None,
                'verbose': False,
                'runall': False,
                'autopick': None,
                'station': False,
                'stream_order': 'same',
                'data_label': None,
                'output': None,
                'fig_size': (12, 8),
                'line_width': 0.825,
                'border_width': 1,
                'tick_width': 1
                }

    fs = 10
    if 'minispec' in args:
        if args['minispec'] == 'pretty':
            fs = 14
    defaults['fs'] = fs

    for key, value in defaults.items():
        if key not in args:
            args[key] = value

    return args


def check_input(data, args):

    if 'minispec' in args:
        if 'station' not in args:
            if type(args['minispec']) == list or args['minispec'] == 'pretty':
                msg = 'List of plots or minispec="pretty" '
                msg += 'only supported when "station" is provided'
                raise IOError(msg)

    args = get_defaults(args)

    # width = args['width']
    cmap = args['cmap']
    show_modes = args['modes']
    syn = args['syn']
    if not isinstance(cmap, type(iter([]))):
        if (
            not hasattr(cm, cmap)
            and cmap != 'BlackGreysRed'
            and cmap != 'BlackRedGreys'
            and cmap != 'GreensBlues'
            and cmap != 'BlueBlackGreysRed'
            and cmap != 'Grays'
            and cmap.lower() != 'black'
             ):
            msg = "\n\033[93mUnknown colormap. Set to 'rainbow'\033[0m\n"
            print(msg)
            cmap = 'rainbow'

    check_data_files = True
    if type(data) == obspy.core.stream.Stream:
        check_data_files = False
    else:
        if type(data) is not list:
            data = [data]
        data_in = data[:]

    if type(syn) == obspy.core.stream.Stream:
        check_data_files = False
    else:
        if syn:
            if type(syn) is not list:
                syn = [syn]
            syn_in = syn[:]
        else:
            syn_in = None

    d_inds = s_inds = None

    if check_data_files is True:
        try:
            while True:
                if data is None:
                    msg = "Enter:\ndata file location ('PATH/FILE.ahx')\n -->  "
                    data = [input(msg)]
                    data_in = data[:]

                if syn_in is not None and syn is None:
                    msg = "Enter one of the following:\n"
                    msg += "- syn file location: PATH/FILE.ahx\n"
                    msg += "- list of syn files:\n"
                    msg += "\t'PATH1/FILE1.ahx', 'PATH2/FILE2.ahx'\n"
                    msg += "- None\n -->  "
                    syn = input(msg)
                    syn = syn.split(',')
                    syn_in = syn[:]
                    if syn in [['None'], ['none']]:
                        syn = None

                data, syn = files_exists([data, syn])

                d_inds = [i for i, x in enumerate(data) if x is False]
                if syn is not None:
                    s_inds = [i for i, x in enumerate(syn) if x is False]
                else:
                    s_inds = None

                if d_inds:
                    msg = "\n\033[91mData file(s) not found:\033[0m\n"
                    for i in d_inds:
                        msg += "%s\n" % data_in[i]
                        data.remove(False)
                    print(msg)
                    if len(data) == 0:
                        data = None
                    continue

                if s_inds:
                    msg = "\n\033[91mSyn file(s) not found:\033[0m\n"
                    for i in s_inds:
                        msg += "%s\n" % syn_in[i]
                        syn.remove(False)
                    print(msg)
                    if len(syn) == 0:
                        syn = None
                    continue

                break
        except KeyboardInterrupt:
            print('interrupted')
            pass
        args['cmt_id'] = data[0].split('/')[-1].split('.')[0]
        args['file_id'] = '.'.join(data[0].split('/')[-1].split('.')[:-1])
    else:
        args['cmt_id'] = args['cmt']
        args['file_id'] = None
    args['syn'] = syn

    if show_modes is True or show_modes == 'all':
        modes = read_modes()
    elif type(show_modes) is list:
        modes_all = read_modes()
        modes = Modes()
        for mode in show_modes:
            modes += modes_all.select(name=mode)[0]
    elif type(show_modes) is Modes:
        modes = show_modes
    else:
        modes = None
    args['modes'] = modes

    if type(data) == obspy.core.stream.Stream:
        return data, args
    else:
        return data[0], args


def files_exists(check_list):
    """
    Takes a list of files as input, returns the list if thew files exists or
    False for each list entry if not
    param check_list: list of strings
    """
    for i, item in enumerate(check_list):
        if item is None or item is False:
            continue
        if type(item) != list:
            if not os.path.exists(item):
                check_list[i] = False
        else:
            for j, subitem in enumerate(item):
                if subitem is False:
                    continue
                elif not os.path.exists(subitem):
                    check_list[i][j] = False
    return check_list
