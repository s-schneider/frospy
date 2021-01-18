#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 15:49:26 2019

@author: talavera
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import matplotlib.markers as mmarkers

import os
import re
from os.path import join
from obspy.core import AttribDict
from collections import OrderedDict
import time
import pandas as pd

from frospy.postprocessing.read import (query_ev_seg_misfits,
                                      query_inversion_summary,
                                      read_cst_files)
from frospy.core.splittingfunc import Set, SplittingFunc
from frospy.core.splittingfunc import loadmodel # python2 branch
#from frospy.splitting.load import loadmodel # python3 master branch
from frospy.core.splittingfunc.plot import sens_kernel
from frospy.core.database.query import get_db_tables, db_query
from frospy.core.modes import format_name
from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes, Mode
from frospy.core.segment import Segment
from frospy.core.segment import read as read_seg
from frospy.core.setup.settings import read as read_setup
from frospy.core.setup.settings import Setup

from frospy.util.read import read_std_cat, read_std_inv
from frospy.util.base import (format_list, update_progress, split_digit_nondigit,
                            is_number, is_mode, uniq_modes, sort_catalog,
                            max_cc_degrees)
from frospy.plot.nmplt import get_iter_colormap
from frospy.util.read import read_st
from frospy.spectrum.app import spectrum
from frospy.core.spectrum.spectrum import Spectrum
from frospy.core.spectrum.plot import plot_modes, plot_segments

import glob
import obspy


def plot_all_data(spectra=None, rundir=None, cmt=None, comp=None,
                  station=None, seg_suffix='dat', select_by='weight',
                  segment_file=None, Nmin=0, Nmax=10, threshold=None,
                  indices=None,xlim=None,
                  verbose=False, weight_fw=None, color=None, ylim=None,
                  data_path=None, removed_segments=None):

    if data_path is None:
        data_path = '/net/home/talavera/eejit/data/VHZ'

    if cmt is not None:
        segment_file = os.path.join(rundir, "%s.%s" % (cmt, seg_suffix))
        if station is None:
            out = spectrum(data_path, seg_reference=True, minispec='amp',
                           segment_file=segment_file,
                           modes='all', runall=True, plot=True)
        elif station == 'search':
            out = spectrum(data_path, seg_reference=True, minispec='amp',
                           segment_file=segment_file,
                           modes='all')
        else:
            out = spectrum(data_path, seg_reference=True, minispec='amp',
                           segment_file=segment_file,
                           modes='all', station=station)

    # Else plot all files
    else:
        if removed_segments is None:
            removed_segments = []
        modes = read_modes()
        fig, ax = plt.subplots()
        if ylim is not None:
            ax.set_ylim(ylim)
        fig.set_size_inches(12, 8)
        l_height = 0.11
        l_width = 0.012

        if verbose:
            print(rundir)
            print(data_path)

        files = glob.glob("%s/???????/???????.%s*" % (rundir, seg_suffix))
        if len(files) == 0:
            files = glob.glob("%s/???????.%s*" % (rundir, seg_suffix))
        if len(files) == 0:
            print('No files found for: %s/*%s*' % (rundir, seg_suffix))
            return

        Nfiles = len(files)
        fw = None

        all_w = []   #
        lines = {}
        all_seg = Segment()
        if spectra is None:
            all_spectra = []
        else:
            spectra_iter = iter(spectra)
            all_spectra = spectra

        for i, seg_f in enumerate(files):
            if verbose is False:
                update_progress((i+1)/float(Nfiles), 'Loading   Data')
            seg = read_seg(seg_f)
            if verbose is True:
                print(seg_f)
                print(seg)

            if len(seg) == 0:
                continue

            if fw is None:
                fw = [seg[0].fw1-0.075, seg[0].fw2+0.075]

            cmt = os.path.basename(seg_f).split('.')[0]

            if comp[0] == 'Z':
                file = os.path.join(data_path, '%s.ahx' % (cmt))
            else:
                file = os.path.join(data_path,
                                    '%scomponent/%s.ahx' % (comp[0], cmt))

            if spectra is None:
                st = read_st(file)
                if verbose is True:
                    print(file)

            for pick in seg:
                if spectra is None:
                    tr = st.select(station=pick.stats.station)
                    if len(tr) == 0:
                        msg = 'pick not in ahx:'
                        msg += ' removing: %s from %s' % (pick, seg_f)
                        print(msg)
                        seg.remove(pick)
                        seg.write(seg_f, overwrite=True)
                        continue
                    else:
                        tr = tr[0]

                    spec = Spectrum(tr, pick.tw1, pick.tw2)
                    all_spectra += [spec]
                else:
                    spec = next(spectra_iter)

                if select_by == 'weight':
                    if weight_fw is not None:
                        all_w.append(spec.weight(weight_fw[0], weight_fw[1]))
                    else:
                        all_w.append(pick.weight)
                elif select_by == 'peak':
                    if weight_fw is not None:
                        ind = [spec.flabel(weight_fw[0]),
                               spec.flabel(weight_fw[1])]
                    else:
                        ind = [spec.flabel(pick.fw1),
                               spec.flabel(pick.fw2)]
                    peak = spec.data.fft.Data[ind[0]:ind[1]]
                    all_w.append(max(abs(peak)))

                if threshold is not None:
                    if all_w[-1] < threshold:
                        if removed_segments is not None:
                            if pick not in removed_segments:
                                spec.plot(fw[0], fw[1], xlabel='f(mHz)', ax=ax)
                                lines[pick] = ax.lines[-1]
                        else:
                            spec.plot(fw[0], fw[1], xlabel='f(mHz)', ax=ax)
                            lines[pick] = ax.lines[-1]
                        
                    else:
                        if pick not in removed_segments:
                            removed_segments.append(pick)
                else:
                    spec.plot(fw[0], fw[1], xlabel='f(mHz)', ax=ax)
                    lines[pick] = ax.lines[-1]
                all_seg += pick

        # Sort segments, with all_w as reference and
        # convert back to list
        _all_seg = []

        for i, w in enumerate(all_w):
            # print(w, all_seg.picks[i].station, all_spectra[i].stats)
            _all_seg.append([w, all_seg.picks[i]])

        _all_seg.sort(key=lambda x: x[0])
        all_seg = Segment(np.array(_all_seg).transpose()[1])
        all_seg = [x for x in reversed(all_seg.picks)]

        all_w = [x for x in reversed(all_w)]

        if removed_segments is not None:
            for p in removed_segments:
                if p in lines:
                    if lines[p] in ax.lines:
                        ax.lines.remove(lines[p])

        print("")
        print("Biggest peaks:")

        if threshold is None:
            if indices is None:
                print('indices')
                for i, pick in enumerate(all_seg[Nmin:Nmax]):
                    print("%i: %s" % (i+Nmin, pick.stats))
                    removed_segments = get_line(lines, pick, ax, color, comp,
                                                removed_segments)
            else:
                print('Nmax')
                if type(indices) is not list:
                    indices = [indices]
                for i in indices:
                    pick = all_seg[i]
                    print("%i: %s" % (i, all_seg[i].stats))
                    removed_segments = get_line(lines, pick, ax, color, comp,
                                                removed_segments)
        else:
            for i, p isn enumerate(removed_segments):
                print("%i: %s" % (i, p.stats))
        if xlim is not None:
            ax.set_xlim(xlim)
        plot_modes(spec, fw[0], fw[1], modes, ax, l_height, l_width)
        plot_segments(spec, seg, fw[0], fw[1], [ax])

    return removed_segments, all_spectra

#from frospy.postprocessing.plot import plot_all_data # python3
#rundir='/net/home/talavera/eejit/splitting/02s00-07s02/allevents/cst+d00+cc/cst_02s00_c20=15'
#rundir="/net/home/talavera/eejit/splitting/03s00-08s02-09s02/allevents/cst+d00+cc/08s02-09s02/03s00"
#rundir="/net/home/talavera/eejit/splitting/04s00-10s02-11s02/goodevents/cst+d00+cc/10s02-11s02/04s00"
#rundir="/net/home/talavera/eejit/splitting/05s00-13s02/deepevents/cst+d00+cc/cst_05s00_13s02/05s00"
#rundir="/net/home/talavera/eejit/splitting/01s00-04s02-00s10/allevents/cst+d00+cc/good/01s00"
rundir="/net/home/talavera/eejit/splitting/11s00-27s02/allevents/cst+d00+cc/cst_27s02+c00=20/"
seg, spectra = plot_all_data(spectra=None, rundir=rundir, 
                             comp='Z',seg_suffix='dat', Nmax=0, color='red',
                             verbose=False,)
#                             xlim=[2.48,2.54],#)
#                             ylim=[0,5e-6])

#replace spectra=None to spectra when run once already