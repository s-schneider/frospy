#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import glob
import os
import re
import sys
import numpy as np

import obspy
try:
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm
    from frospy.plot.nmplt import (format_exponent,
                                 plot_speclat)
    from frospy.util.picker import pick_segment_window
    import cartopy
    import cartopy.crs as ccrs
    NOWINDOW = False
except ImportError:
    NOWINDOW = True
    pass


from frospy.util.base import (pairwise,
                            get_stations_from_stream,
                            find_peaks, get_times, fourier_transform,
                            get_local_extrema, signal2noise, signal2fwhm,
                            takeClosest, group_consecutive_numbers,
                            mask_data, neighbouring_minima)

from frospy.util.read import read_cmt_file, read_std_cat, read_st
from frospy.core.segment import read as read_segments
from frospy.core.segment import Segment
from frospy.core.pick import Pick
from frospy.core.modes import Modes
from frospy.core.modes import read as read_modes
from frospy.core import Spectrum
from frospy.spectrum.help import printhelp

from obspy.core import AttribDict
from obspy import Stream

try:
    from scipy.interpolate import interp1d
    SCIPY = True
except ImportError:
    SCIPY = False

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass


def set_plot_args(main, spec):
    args = {'spec': spec, 'fw': main.fw, 'modes': main.modes,
            'segments': main.segments,
            'st_work': main.st_work,
            'inv': main.inv, 'cmt': main.cmt, 'cmt_id': main.cmt_id,
            'noisewin': main.noisewin, 'cmap': main.cmap,
            'cmap_highlight': main.cmap_highlight,
            'fig': main.rfig, 'ax_old': main.rax, 'seg_ax': main.seg_ax,
            'fs': main.fs, 'zoom_sl': main.zoom_sl, 'width': main.line_width,
            'overlap': main.overlap, 'minispec': main.minispec,
            'save_fig': main.savefig, 'fig_size': main.fig_size}
    plot_args = AttribDict(args)
    return plot_args


def get_misfit_snr(spectrum, pick):

    fwhm = 0  #spectrum.signal2fwhm(pick, max_peak='data')

    if spectrum.syn is None:
        return fwhm, None, None

    # If Fxx is not None continue with the calculation
    rmf = []
    cmf = []
    # misfit is calculated for all synthetics, currently picking only the first
    rmf_t, cmf_t = spectrum.misfit(pick)[0]
    rmf.append(rmf_t)
    cmf.append(cmf_t)

    return fwhm, cmf, rmf


def get_pick(spectrum, fw, modes=None, weighting='sum', event=''):
    weight = spectrum.weight(fw[0], fw[1], rule=weighting)
    tw1 = spectrum.stats.tw[0] / 3600.
    tw2 = spectrum.stats.tw[1] / 3600.

    header = {'station': spectrum.stats.station.code,
              'channel': spectrum.stats.record.channel,
              'event': event}

    pick = Pick(header=header, fw1=fw[0], fw2=fw[1], tw1=tw1, tw2=tw2,
                weight=weight)

    seg_modes = Modes()
    if modes is None:
        modes = read_modes()
    for m in modes:
        if m.freq > fw[0] and m.freq < fw[1]:
            seg_modes += m
    pick.stats.modes = seg_modes
    pick.stats.snr = spectrum.signal2noise(pick)[0]
    return pick


def get_spectrum(tr, st_syn, tw, i, taper_shape, syn_label=None,
                 stream_order='same'):

    station = tr.stats.station
    if st_syn is not None:
        syn_traces = Stream()
        if stream_order == 'same':
            for s in st_syn:
                try:
                    syn_traces += s[i]
                except IndexError:
                    continue
        else:
            for s in st_syn:
                try:
                    syn_traces += s.select(station=station)
                except IndexError:
                    continue
    else:
        syn_traces = None
    spec = Spectrum(tr, tw[0], tw[1], syn_stream=syn_traces,
                    syn_label=syn_label, taper=taper_shape)

    return spec


def quit_spectrum(main, save_files=True):

    if main.file_id:
        if save_files:
            if main.segments:
                segfile = "%s.pickle" % main.cmt_id
                main.segments.write(segfile)
                main.segments.write(main.cmt_id, format='segment')

            if main.save_stream:
                stfile = "%s.ahx" % main.cmt_id
                write_stream(stfile, main.st)

        msg = "\033[92mRemoving tmp files\033[0m"
        print(msg)
        fname = main.file_id + '*' + '.tmp*'
        cwd = os.getcwd()
        cwd = cwd + '/' + fname
        files = glob.glob(cwd)
        for f in files:
            os.remove(f)
    return


def init_time_window_test(main, kwargs, vpoint):
    tr, syn_trs, seg, fw = main.tr, main.syn_trs, main.segments, main.fw
    taper_shape, modedata, tw = main.taper_shape, main.modes, main.tw
    synlabel = main.synlabel
    max_peak = 'data'
    interp = False
    peak_order = 1
    plot_tw = False
    snr_tw = False
    times = None
    mode = None
    verbose = False

    if None in kwargs:
        plot_tw = True

        times = []
        mode = []
        ans = input('Input starttime in h -->  ')
        times.append(float(ans))

        ans = input('Input endtime(s) in h (sep by ",")-->  ')
        if len(ans.split(',')) == 2:
            for a in ans.split(','):
                times.append(float(a))
        else:
            msg = 'Wrong format, new input:\n'
            print(msg)
            ans = input('Input endtime 1 in h -->  ')
            times.append(float(ans))

            ans = input('Input endtime 2 in h -->  ')
            times.append(float(ans))

        ans = input('Input stepsize in h -->  ')
        times.append(float(ans))

        msg = 'Weighting of peaks (Integer: 1, 2, ...) -->  '
        ans = input(msg)
        peak_order = float(ans)

        msg = 'Test for certain mode? (y/n)-->  '
        ans = input(msg)

        if ans in ['y', 'Y', 'yes', 'Yes']:
            msg = 'Enter Mode name -->  '
            ans = input(msg)
            mode = get_mode(ans, modedata)

            msg = 'Create envelope based on peaks in Amp.-Spect?'
            msg += ' (y,n) -->  '
            ans = input(msg)
            if ans in ['y', 'Y', 'yes', 'Yes']:
                interp = True

                msg = 'maximum peak based on interpolated data? '
                msg += '(y,n) -->  '
                ans = input(msg)
                if ans in ['y', 'Y', 'yes', 'Yes']:
                    max_peak = 'interpolated'

            msg = 'Calculate SNR * SNFWHMR?'
            msg += ' (y,n) -->  '
            ans = input(msg)
            if ans in ['y', 'Y', 'yes', 'Yes']:
                snr_tw = 'true'

        ans = input('Save figure? (y/n) -->  ')
        if ans in ['y', 'Y']:
            save = True
        else:
            save = False

    elif len(kwargs) >= 4:
        try:
            times = [float(kwargs[0]), float(kwargs[1]),
                     float(kwargs[2]), float(kwargs[3])]
            save = False
            for kw in kwargs[4:]:
                if kw == 'save':
                    save = True
                elif kw.isdigit():
                    peak_order = int(kw)
                elif kw == 'interp':
                    interp = True
                elif kw == 'interpolatedpeak':
                    max_peak = 'interpolated'
                elif kw == 'plot':
                    plot_tw = True
                elif kw == 'snr':
                    snr_tw = True
                elif kw == 'verbose':
                    verbose = True
                else:
                    mode = get_mode(kw, modedata)
        except Exception:
            msg = "Invalid input. \n"
            msg += "Error message: %s" % sys.exc_info()[1]
            print(msg)

    else:
        msg = 'Usage: \n'
        msg += 'tw test start - t_start_1 t_start_2'
        msg += ' t_end t_step (kwargs)\n'
        msg += 'or\n'
        msg += 'tw test end -  t_start t_end_1 t_end_2 t_step (kwargs)\n'
        print(msg)
        return

    # try:
        if verbose:
            msg = '\nPerforming timewindow test with following parameters:\n'
            msg += (
                'max_peak = %s\
                \ninterp = %s\
                \npeak_order = %s\
                \nplot_tw = %s\
                \nsnr_tw = %s\
                \ntimes = %s\
                \nmode = %s\n'
                )
            msg = msg % (max_peak, interp, peak_order, plot_tw, snr_tw,
                         times, mode)
            print(msg)

    dtitle = "%s %s" % (tr.stats.station, tr.stats.channel)
    # try:
    tw = time_window_test(tr, taper_shape, fw, times,
                          modedata, mode, title=dtitle,
                          savefig=save, vpoint=vpoint,
                          max_peak=max_peak, interp=interp,
                          peak_order=peak_order, plot=plot_tw,
                          snr_tw=snr_tw, segment=seg,
                          weighting='integrate')

    if synlabel:
        stitle = synlabel[0]
        ts = syn_trs[0]
        time_window_test(ts, taper_shape, fw, times, modedata,
                         mode, title=stitle, savefig=save,
                         vpoint=vpoint, max_peak=max_peak, interp=interp,
                         peak_order=peak_order, plot=plot_tw,
                         snr_tw=snr_tw,
                         segment=seg, weighting='integrate')

    # except Exception:
    #         msg = "Couldn't perform timewindow test. \n"
    #         msg += "Error message: %s" % sys.exc_info()[1]
    #         print(msg)
    #         return tw

    if snr_tw:
        return [tw[0] * 3600., tw[1] * 3600.]
    else:
        return tw


def call_taper_test(tr, trsyn, wstart, wend, tw, modedata, kwargs):
    if None in kwargs:
        savefig = False
    else:
        try:
            savefig = kwargs[0]
        except IndexError:
            savefig = False

    taper_test(tr, trsyn, wstart, wend, tw, modedata, savefig)
    return


def taper_test(tr, trsyn, wstart, wend, tw, modedata, savefig):

    title = 'Taper shape test'
    taper_list = ['hanning', 'bartlett', 'hamming', 'blackman', 'kaiser',
                  'boxcar']

    colormap = get_iter_colormap(taper_list, 'rainbow')

    fig, ax = plt.subplots()

    for shape in taper_list:
        c = next(colormap)
        times = get_times(tr, tw[0], tw[1], None, None, verbose=False)
        tstart, nstart, nend, t, Htstart, Htend, Terr = times

        f, Fxx, delomeg = fourier_transform(tr, Htstart, Htend,
                                            nstart, nend, shape)

        startlabel = get_freq_index(wstart, delomeg)
        endlabel = get_freq_index(wend, delomeg)

        label = '%s' % shape
        ax.plot(f[startlabel:endlabel+1], abs(Fxx[startlabel:endlabel+1]),
                label=label, linestyle='dashed', color=c)

    f_lim = [f[startlabel], f[endlabel+1]]
    _plot_modes_in_ax(ax, modedata, f_lim)
    ax.set_xlim(f[startlabel], f[endlabel+1])
    ax = format_exponent(ax)
    ax.legend()
    plt.suptitle(title)

    if savefig:
        fig.set_size_inches(12, 8)
        fig.savefig('%s_taper_test.png' % (tr.stats.station), dpi=300,
                    orientation='landscape')
    plt.show()
    return


def time_window_test(tr, taper_shape, fw, times, modedata, mode,
                     title=None, savefig=False,
                     vpoint='end', max_peak='data', interp='True',
                     peak_order=2., plot=True, snr_tw=False, segment=None,
                     weighting='integrate'):

    """
    Documentation follows
    """
    if segment:
        seg = segment.copy()
    else:
        seg = None

    wstart = fw[0]
    wend = fw[1]
    # Getting time range from time values stored in mode
    window = np.array([times[0], times[1], times[2], times[3]]) * 3600.
    gtimes = get_times(tr, window[0], window[2], None, None, verbose=False)
    tstart, nstart, nend, t, Htstart, Htend, Terr = gtimes
    tend = np.round(Htend)
    tq = None

    if vpoint == 'end':
        test_times = np.arange(window[1], tend, window[3])
        cnorm = plt.Normalize(vmin=int(window[1] / 3600.),
                              vmax=int((tend) / 3600.))
    else:
        test_times = np.arange(window[0], window[1], window[3])
        cnorm = plt.Normalize(vmin=int(window[0] / 3600.),
                              vmax=int((window[1]+window[3]) / 3600.))
    s2fwhm = []
    snr = []

    t_old = 0
    final_times = []
    Farr = []

    if plot:
        fig, ax = plt.subplots()
        ax = range(2)
        ax[0] = plt.subplot2grid((2, 1), (0, 0))
        ax[1] = plt.subplot2grid((2, 1), (1, 0))
        cax = plt.cm.ScalarMappable(cmap='rainbow', norm=cnorm)
        cax._A = []
        colormap = get_iter_colormap(test_times, 'rainbow')
        clim = []

    for _t in test_times:
        if plot:
            c = next(colormap)
        if vpoint == 'end':
            gtimes = get_times(tr, window[0], _t, None, None, verbose=False)
        else:
            gtimes = get_times(tr, _t, window[2], None, None, verbose=False)

        tstart, nstart, nend, t, Htstart, Htend, Terr = gtimes

        if Htend == t_old and vpoint == 'end':
            continue
        f, Fxx, delomeg = fourier_transform(tr, Htstart, Htend,
                                            nstart, nend, taper_shape)

        startlabel = get_freq_index(wstart, delomeg)
        endlabel = get_freq_index(wend, delomeg)

        # calculate signal to full widht half magnitude ratio
        if mode is not None:
            if len(mode) == 1:

                seg = auto_set_pick(tr, f, Fxx, delomeg, modedata,
                                    Htstart, Htend, None, None, seg,
                                    taper_shape, weighting, mode,
                                    verbose=False, single_seg=True)
                pick = seg.select(station=tr.stats.station)[0]
                Q = mode.Q
                f_center = mode.freq
                Hqc = 1.1 * Q / (f_center/1000.) / 3600.
                i_center = int(f_center * 2. *
                               np.pi/(1000. * delomeg) + 1
                               - startlabel
                               )

                seg_ind1 = int(pick.fw1*2.*np.pi/(1000.*delomeg))+1
                seg_ind2 = int(pick.fw2*2.*np.pi/(1000.*delomeg))+1

                s2f = signal2fwhm(Fxx[seg_ind1-1:seg_ind2+2],
                                  i_center, max_peak=max_peak,
                                  interpolate=interp,
                                  peak_order=peak_order)
                s2fwhm.append([[Htstart/3600., Htend/3600.], s2f])

                if snr_tw:
                    snr_out = signal2noise(Fxx, pick, delomeg)
                    snr_t, nwins, ol_err = snr_out
                    snr.append([[Htstart/3600., Htend/3600.], snr_t])
        else:
            Hqc = None

        if plot:
            if np.round(_t) == Hqc:
                if vpoint == 'end':
                    l_t = (_t - Htstart)/3600.
                else:
                    l_t = (_t)/3600.
                label = r'%.1f h (qcycle)' % (l_t)
                ax[0].plot(f[startlabel:endlabel+1],
                           abs(Fxx[startlabel:endlabel+1]),
                           label=label, linestyle='solid', color='black')
            else:
                if vpoint == 'end':
                    l_t = (_t - Htstart)/3600.
                else:
                    l_t = (_t)/3600.
                label = r'%.1f h' % (l_t)
                ax[0].plot(f[startlabel:endlabel+1],
                           abs(Fxx[startlabel:endlabel+1]),
                           label=label, linestyle='dashed', color=c)
            clim.append(l_t)

        Farr.append(abs(Fxx[startlabel:endlabel+1]))
        final_times.append(_t)
        t_old = Htend

    s2fwhm = np.array(s2fwhm)
    if s2fwhm.any():
        if snr_tw:
            snr = np.array(snr)
            sms = snr.transpose()[1] * s2fwhm.transpose()[1]
            tw_pick = s2fwhm[sms.argmax()][0]
            max_value = s2fwhm[sms.argmax()][1]
        else:
            tw_pick = s2fwhm[s2fwhm.transpose()[1].argmax()][0]
            max_value = s2fwhm[s2fwhm.transpose()[1].argmax()][1]
        f_pick = int(f_center * 2. * np.pi/(1000. * delomeg)) + 1

        msg = "\n\tSuggested timewindow:\t%.1f - %.1f hrs" % (tw_pick[0],
                                                              tw_pick[1])
        msg += "\n\tPeak Value:\t%.1e\n" % max_value
        qc_fac = tw_pick[1]/Hqc
        msg += "\tEndtime is:\t%.2f x Qcycle (%.2f hrs)\n" % (qc_fac, Hqc)
        print(msg)
    else:
        tw_pick = [None, None]

    if plot:
        cbar = fig.colorbar(cax, ax=ax[0], orientation="vertical")

        f_lim = [f[startlabel], f[endlabel+1]]
        _plot_modes_in_ax(ax[0], modedata, f_lim)
        _plot_modes_in_ax(ax[1], modedata, f_lim, 'white')

        if vpoint == 'end':
            if s2fwhm.any():
                tw_picked = tw_pick[1]
            cbarlabel = r"$t_{end}$ in h"
            extent = [f[startlabel], f[endlabel+1],
                      int((final_times[0] - Htstart) / 3600.),
                      int((final_times[-1] - Htstart) / 3600.)]

        elif vpoint == 'start':
            if s2fwhm.any():
                tw_picked = tw_pick[0]
            cbarlabel = r"$t_{start}$ in h"
            extent = [f[startlabel], f[endlabel+1],
                      int((final_times[0])/3600.),
                      int((final_times[-1])/3600.)]

        cbar.set_label(cbarlabel)

        ax[0].set_xlim(f[startlabel], f[endlabel+1])

        im = ax[1].imshow(Farr, extent=extent, origin='lower', aspect='auto')

        if tq:
            ax[1].annotate('qcycle', fontsize=10, xy=(f[startlabel+1],
                           (tq - Htstart)/3600.),
                           xycoords='data', xytext=(0, 0),
                           textcoords='offset points',
                           color='white',
                           )
        if s2fwhm.any():
            ax[1].plot(f[f_pick], tw_picked, 'v', color='white')

        cbarim = fig.colorbar(im, ax=ax[1], orientation="vertical")
        cbarim.set_label(r"Amplitude")
        ax[1].set_xlim(f[startlabel], f[endlabel+1])
        ax[1].set_ylabel(cbarlabel)
        ax[1].set_xlabel('frequency (mHz)')

        if vpoint == 'end':
            plt.suptitle(r"%s  $t_{start}$: %i h fixed" % (title,
                                                           Htstart/3600.))
        else:
            plt.suptitle(r"%s  $t_{end}$: %i h fixed" % (title, Htend/3600.))
        if savefig:
            fig.set_size_inches(12, 8)
            fig.savefig('%s-%s_tw_test.png' % (tr.stats.station, title),
                        dpi=300, orientation='landscape')
        else:
            plt.show()

    return [tw_pick[0], tw_pick[1]]


def report_input(main):
    msg = '\n'
    if main.st and main.st_syn and main.cmt and main.segments:
        msg += "\t All files succesfully read\n"
    if not main.st:
        msg += "\t No data files loaded\n"
    if not main.st_syn:
        msg = "\t No synthetic files loaded\n"
    if not main.cmt:
        msg += "\t No cmt file loaded\n"
    if not main.segments:
        msg += "\t No segment file read.\n"

    return msg


def list_traces(st):
    for i, tr in enumerate(st):
        print("%i: %s" % (i, tr))
    return


def load_segments(kwargs, verbose=True):
    if type(kwargs) is not list:
        return None
    elif type(kwargs[0]) is Segment:
        return kwargs[0]
    elif len(kwargs) > 1 or None in kwargs:
        file = input('Enter filename\n -->  ')
    elif len(kwargs) == 1:
        file = kwargs[0]

    try:
        segments = read_segments(file)
        msg = "\033[92m\t {Nseg} segment(s) loaded\033[0m"
        msg = msg.format(Nseg=len(segments))
        if verbose:
            print(msg)
        return segments
    except IOError:
        msg = "\t segment-file not found!"
        print(msg)
        return None
    except ValueError:
        msg = "\t segment-file in wrong format, must be: \n\n"
        msg += "Station Name 1\n"
        msg += "freq1 freq2\n"
        msg += "time1 time2\n"
        msg += "weighting\n"
        msg += "Station Name 2\n"
        msg += "freq1 freq2\n"
        msg += "time1 ...\n"
        print(msg)
        return None
    except Exception:
        msg = "Unexpected error! \n"
        msg += "Error message: %s: %s" % (sys.exc_info()[0],
                                          sys.exc_info()[1])
        print(msg)
        return None


def load_cmt(id=None, file=None, verbose=True):
    # msg = "Loading cmt ..."
    # print(msg)
    if file:
        try:
            cmt = read_cmt_file(file)
            if verbose:
                print("\t cmt loaded")
        except IOError:
            if verbose:
                print("\t cmt not found, continue")
            cmt = None
    else:
        try:
            cat = read_std_cat(id)
            cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
            m0 = cat[0].focal_mechanisms[0].moment_tensor.scalar_moment
            mw = cat[0].magnitudes[0].mag
            depth = cat[0].origins[0].depth / 1000.
            cmt = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp,
                   mw, m0, depth]
            if verbose:
                print("\t cmt loaded")
        except IOError:
            if verbose:
                print("\t cmt not found, continue")
            cmt = None
        except TypeError:
            if verbose:
                print("\t cmt not in catalog, continue")
            cmt = None
    return cmt


def read_data(data, main, verbose=True):
    if type(data) == obspy.core.stream.Stream:
        st = data
    else:
        st = read_st(data, 'ah')
    # st_tmp.sort(['station'])
    if main.verbose:
        msg = "\033[92m\t {Ndata} trace(s) loaded\033[0m".format(Ndata=len(st))
        print(msg)

    if main.segments is not None:
        st_work = Stream()
        for pick in main.segments:
            try:
                st_work += st.select(station=pick.stats.station)[0]
            except IndexError:
                continue
    else:
        st_work = st.copy()

    main.st = st
    main.st_work = st_work
    main = read_syn(main, verbose)

    return main


def read_syn(main, verbose=True):
    syn, st, segments = main.syn, main.st, main.segments
    stat_d = get_stations_from_stream(st)
    st_syn = []
    st_syn_work = []
    if type(syn) is str or type(syn) == obspy.core.stream.Stream:
        if type(syn) == obspy.core.stream.Stream:
            stsyn_tmp = syn
        else:
            stsyn_tmp = read_st(syn)

        stat_s = get_stations_from_stream(stsyn_tmp)
        # check if at least one matchin station is in syn
        if len(set(stat_s).intersection(stat_d)) < 1:
            msg = "No matching stations between data and syntheics in:\n"
            msg += "%s" % file
            raise IOError(msg)

        if segments is not None:
            st = Stream()
            for pick in segments:
                st += stsyn_tmp.select(station=pick.stats.station)[0]

        else:
            st = stsyn_tmp

        #stsyn_tmp.sort(['station'])
        st_syn.append(stsyn_tmp)
        st_syn_work.append(st)

        if main.verbose:
            "\033[92m\t {Ndata} trace(s) loaded\033[0m"
            msg = "\033[92m\t ({N0}/{N}) synthetics loaded\033[0m"
            print(msg.format(N0=1, N=len(st_syn)))

    elif type(syn) is list:
        i = 0
        for file in syn:
            stsyn_tmp = read_st(file)
            stat_s = get_stations_from_stream(stsyn_tmp)
            # check if at least one matchin station is in syn
            if len(set(stat_s).intersection(stat_d)) < 1:
                msg = "No matching stations between data and syntheics in:\n"
                msg += "%s" % file
                print(msg)
                continue
            else:
                if segments is not None:
                    st = Stream()
                    for pick in segments:
                        st += stsyn_tmp.select(station=pick.stats.station)[0]
                else:
                    st = stsyn_tmp
                i += 1
                #stsyn_tmp.sort(['station'])
                st_syn.append(stsyn_tmp)
                st_syn_work.append(st)
                if verbose:
                    print("\t (%i/%i) synthetics loaded" % (i, len(syn)))

    elif syn is None:
        st_syn, st_syn_work = None, None

    else:
        msg = "Wrong input format of synthetic data files."
        msg += "Must be string or list of strings"
        raise IOError(msg)

    main.st_syn, main.st_syn_work = st_syn, st_syn_work
    return main


def show_picks():
    uans = input('Show picks? (y/n) \n')
    if uans == 'y':
        ans = True
    else:
        ans = False
    return ans


def append_segment(spec, seg_picks, segments=None, weighting='integrate',
                   cmt=None):

    allmodes = read_modes()
    for sp in pairwise(seg_picks):
        seg_modes = Modes()
        for m in allmodes:
            if m.freq > sp[0] and m.freq < sp[1]:
                seg_modes += m
        pick = Pick(header={'station': spec.stats.station.code,
                            'channel': spec.stats.record.channel,
                            'event': cmt})
        pick.fw1 = sp[0]
        pick.fw2 = sp[1]
        pick.tw1 = spec.stats.tw[0] / 3600.
        pick.tw2 = spec.stats.tw[1] / 3600
        pick.weight = spec.weight(sp[0], sp[1])
        pick.stats.modes = seg_modes
        if not segments:
            segments = Segment(pick)
        else:
            segments += pick

    return segments


def write_stream(fname, st, overwrite=False):
    try:
        if not overwrite:
            if os.path.exists(fname):
                fname += '.new'
                msg = '\033[93mStream-file exist, not overwriting\033[0m'
                print(msg)
        filename = fname
        st.write(filename, format='AH')
        msg = "\033[92mStream-file written to %s\033[0m" % filename
        print(msg)
    except IOError:
        msg = "\033[91mCan't save file\n"
        msg += "Error message: %s\033[0m" % sys.exc_info()[1]
        print(msg)
    return


def print_station(main):
    file_id, stream = main.file_id, main.st_work
    wfile = file_id + '.stations'
    stations = get_stations_from_stream(stream)
    with open(wfile, 'w') as fh:
        for station in stations:
            fh.write("%s\n" % station)
    msg = "\033[92mStation file written to %s\033[0m" % wfile
    print(msg)
    return


def get_show_modes(modes, kwargs=None):
    if None in kwargs:
        msg = 'Showing modes? (all / none / names (separated by space)'
        msg += ', e.g. 0S15) -->  '
        ans = input(msg)
    else:
        ans = ' '.join(kwargs)
    if ans == 'y' or ans == 'all':
        return modes

    elif ans == 'n' or ans == 'none' or ans == 'off':
        return None

    else:
        while True:
            try:
                return get_mode(ans, modes)
            except KeyError:
                msg = 'Wrong name of mode, enter mode or "quit" -->  '
                ans = input(msg)
                if ans == 'quit':
                    return None


def get_mode(names, modes):
    if names is None:
        names = input('Which mode? (name, e.g. 0S15) -->  ')
    if type(names) is str:
        names = names.split()

    modes_new = Modes()
    for name in names:
        # Checking if only S or T should be displayed
        if name in ['T', 't', 'S', 's']:
            modes_new = modes.select(mtype=name.upper())
            return modes_new
        # Checking for right format: Digit Letter Digit
        N = re.split('(\d+)', name)[1].lstrip('0')
        M = re.split('(\d+)', name)[2].capitalize()
        L = re.split('(\d+)', name)[3].lstrip('0')
        if len(N) == 0:
            N = '0'
        if len(L) == 0:
            L = '0'

        mode = [N, M, L]
        mode = str().join(mode)
        modes_new += modes.select(name=mode)

    if len(modes_new) == 0:
        modes_new = None

    return modes_new


def search_stream(main, search):
    stream, syn_list, segments = main.st, main.st_syn, main.segments
    st_search = Stream()
    if search[0] == 'segments':
        if segments:
            stations = get_stations_from_stream(stream)
            for p in segments:
                if (
                    p.stats.station in stations and
                    len(st_search.select(station=p.stats.station)) == 0
                   ):
                    st_search.append(stream.select(station=p.stats.station)[0])
    else:
        if type(search) != list:
            search = [search]
        for item in search:
            for trace in stream:
                for attr in trace.stats:
                    if trace.stats[attr] == item:
                        st_search.append(trace)

    if len(st_search) != 0:
        sts_search = []
        if syn_list:
            for sts in syn_list:
                sts_hit = Stream()
                for trace in st_search:
                    sts_hit.append(sts.select(station=trace.stats.station)[0])
                sts_search.append(sts_hit)

        main.st_work = st_search
        main.st_syn_work = sts_search
        main.i_old = main.i
        main.i = 0
        return main

    else:
        main.st_work = main.st
        main.st_syn_work = main.st_syn
        main.i_old = None
        print('Search-term not found')
        return main


def select_channel(stream, streamsyn, chan):
    if len(chan) == 1:
        chan = "*" + chan + "*"

    st_select = stream.select(channel=chan)

    if streamsyn:
        stsyn_select = []
        for st in streamsyn:
            stsyn_select.append(st.select(channel=chan))

        if len(st_select) != len(stsyn_select):
            chan = chan + 's'
            stsyn_select = []
            for st in streamsyn:
                stsyn_select.append(st.select(channel=chan))
    else:
        stsyn_select = None

    return st_select, stsyn_select


def timewindow(tr, kwargs=None):
    check_input = True

    if None not in kwargs:
        try:
            start = float(kwargs[0])
            end = float(kwargs[1])
            check_input = False
        except Exception:
            msg = 'Usage: tw - start end'
            print(msg)
            check_input = True
    else:
        # Htart must be in seconds, input is in Hours, du to practicality
        while True:
            try:
                start = float(input('Input start time in hours  > '))
                end = float(input('Input end time in hours    > '))
                check_input = False
                break
            except KeyboardInterrupt:
                break
            except Exception:
                print('Wrong input')
                continue

    while True:
        if check_input:
            try:
                start = float(input('Input start time in hours  > '))
                end = float(input('Input end time in hours    > '))
            except KeyboardInterrupt:
                break
            except Exception:
                print('Wrong input')
                continue

        if start >= end:
            msg = 'Start must be smaller then end'
            print(msg)
            check_input = True
            continue

        elif type(start) != int and type(start) != float:
            msg = 'Wrong input, must be number'
            print(msg)
            check_input = True
            continue

        elif type(end) != int and type(end) != float:
            msg = 'Wrong input, must be number'
            print(msg)
            check_input = True
            continue
        else:
            check_input = False
            break

    maxtw = tr.stats.endtime - tr.stats.starttime

    end_org = end
    if end > maxtw:
        end = maxtw - 1.
        print('\033[93mEndtime too big, set to %i h\033[0m' % int(end))

    tw = [start, end]
    tw_org = [start, end_org]
    return tw, tw_org


def freqwindow(kwargs=None):
    check_input = False
    fw = [0, 0]
    if None not in kwargs:
        try:
            fw[0] = float(kwargs[0])
            fw[1] = float(kwargs[1])
        except Exception:
            msg = 'Usage: fw - start end'
            print(msg)
            check_input = True
    else:
        while True:
            try:
                fw[0] = float(input('Input start frequency in mHz > '))
                fw[1] = float(input('Input end frequency in mHz   > '))
                break
            except KeyboardInterrupt:
                break
            except Exception:
                print('Wrong input')

    while True:
        if check_input:
            try:
                fw[0] = float(input('Input start frequency in mHz > '))
                fw[1] = float(input('Input end frequency in mHz   > '))
            except KeyboardInterrupt:
                break
            except Exception:
                print('Wrong input')
                msg = 'Wrong input, must be number'

        if fw[0] >= fw[1]:
            msg = 'Start must be smaller then end'
            print(msg)
            check_input = True

        elif type(fw[0]) != int and type(fw[0]) != float:
            msg = 'Wrong input, must be number'
            print(msg)
            check_input = True

        elif type(fw[1]) != int and type(fw[1]) != float:
            msg = 'Wrong input, must be number'
            print(msg)
            check_input = True

        else:
            break

    return fw


def calc_weights(tr, tw, seg_picks, taper_shape, weighting='integrate'):

    freq, amp = get_freq_amp(tr, tw, seg_picks, taper_shape)

    if freq is None and amp is None:
        return None
    else:
        if weighting == 'integrate':
            weight = np.trapz(amp)
        return weight


def get_freq_amp(tr, tw, seg_picks, taper_shape, force_write=True):

    times = get_times(tr, tw[0], tw[1], None, None, verbose=False)
    tstart, nstart, nend, t, tw[0], tw[1], Terr = times

    if not force_write:
        if Terr:
            return None, None

    freq, Fxx, delomeg = fourier_transform(tr, tw, nstart, nend, taper_shape)

    if type(seg_picks) == Pick:
        wstart = seg_picks.fw1
        wend = seg_picks.fw2
    else:
        wstart = seg_picks[0]
        wend = seg_picks[1]

    startlabel = get_freq_index(wstart, delomeg)
    endlabel = get_freq_index(wend, delomeg)

    amp = np.sqrt(Fxx[startlabel:endlabel+1].real**2. +
                  Fxx[startlabel:endlabel+1].imag**2.)

    return freq, amp


def prep_fpeaks(main):

    if not main.segments:
        msg = 'No segments picked, or segment file loaded'
        print(msg)
        return

    if not main.stsyn_work:
        msg = 'No synthetic data loaded'
        print(msg)
        return
    print_fpeaks(main.file_id, main.st_work, main.stsyn_work,
                 main.tw, main.fw, main.cmt, main.segments, main.taper_shape)
    return


def print_fpeaks(file_id, data, syn, tw, fw, cmt, segments, taper_shape):

    if not segments:
        return
    wfile = file_id + '.peaks'

    with open(wfile, 'w') as fh:
        if cmt:
            h1 = "# CMT\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(cmt[0]), str(cmt[1]),
                                                      str(cmt[2]), str(cmt[3]),
                                                      str(cmt[4]), str(cmt[5]))
            fh.write(h1)
        h2 = "# Frequencies in mHz\n"
        h3 = "# station\tfstart\tfend\tf_max_d\ta_max_d\tf_max_s\ta_max_s\n\n"

        fh.write(h2)
        fh.write(h3)

        st_total = len(segments)
        k = 1
        for pick in segments:
            station = pick.stats.station
            trace = data.select(station=station)[0]
            syn_stream = syn[0].select(station=station)
            spec = Spectrum(trace, tw[0], tw[1], syn_stream)
            msg = 'Calculating freq peaks of station '
            msg = msg + '%5s, %i/%i'

            print('', end='\r')
            sys.stdout.flush()
            print(msg % (pick.stats.station, k, st_total), end='\r')
            sys.stdout.flush()

            s = spec.flabel(pick.fw1)
            e = spec.flabel(pick.fw2)

            f = spec.stats.freq[s:e+1][np.abs(
                            spec.data.fft.Data[s:e+1]).argmax()]
            amp = np.abs(spec.data.fft.Data[s:e+1]).max()

            fsyn = spec.stats.freq[s:e+1][np.abs(
                    list(spec.syn[0].fft.values())[0][s:e+1]).argmax()]
            ampsyn = np.abs(
                    list(spec.syn[0].fft.values())[0][s:e+1]).max()

            rds = amp / ampsyn
            entries = "%s\t%f\t%f\t%f\t%E\t%f\t%E\t%E\n"
            fh.write(entries % (pick.stats.station, pick.fw1, pick.fw2, f, amp,
                                fsyn, ampsyn, rds))
            k += 1
    print('\n')
    print('peaks saved to file %s \n' % wfile)
    return


def printw(file_id, f, xp, xpsyn, Fxx, Fxxsyn, stats, tstart, tend, fstart,
           fend, rmf, cmf):

    print('Current data window of station %s saved \n\n' %
          stats.station)
    wfile = file_id + stats.channel + '.window'
    ampfile = file_id + stats.channel + '.af'
    synampfile = file_id + stats.channel + '.afsyn'
    phasefile = file_id + stats.channel + '.pf'
    synphasefile = file_id + stats.channel + '.pfsyn'

    with open(wfile, 'w') as fh:
        fh.write("%s\n" % stats.station)
        if rmf and cmf:
            for r, c in zip(rmf, cmf):
                fh.write("%f\t%f\n" % (r, c))
        else:
            fh.write("None\tNone\n")
        fh.write("%f\t%f\n" % (fstart, fend))
        fh.write("%.1f\t%.1f\n" % (tstart/3600., tend/3600.))

    afile = open(ampfile, 'w')
    safile = open(synampfile, 'w')
    pfile = open(phasefile, 'w')
    spfile = open(synphasefile, 'w')

    for freq, amp, samp, phase, sphase in zip(f, Fxx, Fxxsyn, xp, xpsyn):
            afile.write("%f\t%.6E\t%.6E\n" % (freq, amp.real, amp.imag))
            safile.write("%f\t%.6E\t%.6E\n" % (freq, samp.real, samp.imag))
            pfile.write("%f\t%.6E\n" % (freq, phase))
            spfile.write("%f\t%.6E\n" % (freq, sphase))

    afile.close()
    safile.close()
    pfile.close()
    spfile.close()

    return


def process_input(pick, main, msg=None):
    """
    spec, i, j, tr, st, syn_trs, st_syn, st_work,
    st_syn_work, recalc, startlabel, endlabel,
    cmt_id, rmf, cmf,
    segments, cmt, tw_org,
    tw, fw,
    mfac, i_old, modes, modes_all,
    qwindow, file_id, taper_shape, weighting,
    update_maps, syn_label, pick, inv,
    savefig, fs, min_snr, rfig, rax, segment_file,
    save_stream,
    zoom_sl, cmap):
    """
    main.update_maps = True

    if msg is None:
        msg = 'Waiting for user input (type help for all commands)\n -->  '
        uinput = input(msg).split()
    else:
        uinput = msg.split()

    # Check for flag '-' for keyword arguments
    if len(uinput) == 0:
        uin = ''
        kwargs = [None]
    elif '-' in uinput:
        kwargs = []
        uin = ' '.join(uinput).split('-')[0].strip()
        kwargs_tmp = ' '.join(uinput).split('-')[1:]
        for k in kwargs_tmp:
            kwargs.append(k.strip().split()[:])
        if len(kwargs) == 1:
            kwargs = kwargs[0]
    else:
        uin = ' '.join(uinput)
        kwargs = [None]

    # Checking user input arguments
    if uin == 'next' or uin == 'nsp' or uin == '':
        main.i += 1
        main.j = 0

    elif uin == 'nseg':
        main.j += 1

    elif uin == 'pseg':
        if main.j != 0:
            main.j -= 1
        else:
            main.j = 0
            print('Already at beginning of segments.\n')

    elif uin == 'quit' or uin == 'exit' or uin == 'q':
        if main.segments:
            main.segments = update_segments(['tw'], main.spec, main)
        quit_spectrum(main, main.output)
        return None

    elif uin == 'search':
        if None in kwargs:
            msg = 'Enter station attribute(s) (Attrib1, Attrib2):\n'
            kwargs = input(msg).split(',')

        main = search_stream(main, kwargs)

    elif uin == 'bts':
        main.i = 0

    elif uin == 'psp' or uin == '..':
        if main.i != 0:
            main.i -= 1
            main.j = 0
        else:
            main.i = 0
            main.j = 0
            print('Already at the beginning.\n')

    elif uin == 'fw' or uin == 'Fw':
        main.fw = freqwindow(kwargs)
        main.update_maps = False

    elif uin == 'tw' or uin == 'Tw':
        main.tw, main.tw_org = timewindow(main.tr, kwargs)
        if pick:
            try:
                main.segments = update_segments(['tw'], main.spec, main)
            except Exception:
                pass
        main.qwindow = None
        main.update_maps = False

    elif uin == 'list' or uin == 'l' or uin == 'ls':
        list_traces(main.st_work)
        main.recalc = False

    elif uin == 'list all' or uin == 'la':
        for tr in main.st_work:
            print(tr.stats)

    elif uin == 'list picks' or uin == 'list pick':
        if main.segments:
            for pick in main.segments:
                print(pick)

    elif uin == 'rm dublicate picks' or uin == 'rdp':
        if main.segments:
            main.segments.del_duplicates()

    elif uin == 'printw':
        freq = main.spec.stats.freq
        phase = np.angle(main.spec.data.fft)
        phasesyn = np.angle(main.spec.syn[0].fft)
        amp = abs(main.spec.data.fft)
        ampsyn = abs(main.spec.syn[0].fft)
        printw(main.file_id, freq, phase, phasesyn, amp, ampsyn, main.tr.stats,
               main.tw[0], main.tw[1], main.fw[0], main.fw[1], main.rmf,
               main.cmf)
        main.recalc = False

    elif uin == 'help':
        printhelp()
        main.recalc = False

    elif uin == 'modes':
        try:
            main.modes = get_show_modes(main.modes_all, kwargs)
            if main.modes and len(main.modes) == 1:
                main.fw[0] = main.modes[0].freq - 0.15
                main.fw[1] = main.modes[0].freq + 0.15
            main.update_maps = False

        except Exception:
            msg = 'Invalid mode name!'
            print(msg)

    elif uin == 'mfac':
        msg = 'Enter multiplication factor for synthetics  > '
        main.mfac = float(input(msg))
        main.update_maps = False

    elif uin == 'goto':
        main.i = goto(main.st_work, main.i, kwargs)

    elif uin == 'd' or uin == 'delete' or uin == 'rm':
        if kwargs[0] == 'unfit':
            for tr in main.st:
                times = get_times(tr, main.tw[0], main.tw[1], main.tw_org[0],
                                  main.tw_org[1])
                Terr = times[7]
                if Terr:
                    remove_station(main, save_tmp=False)

        elif kwargs[0] == 'w':
            remove_station(main, save_tmp=False)

        elif kwargs[0] == 'picks' or kwargs[0] == 'pick':
            remove_pick(main.spec, main.segments, kwargs)

        else:
            remove_station(main, save_tmp=False)

    elif uin == 'dall' or uin == 'rmall':
        remove_station(main, save_tmp=False, all=True)

    elif uin == 'save':
        save_data(main, main.file_id, overwrite=True)
        main.recalc = False

    elif uin == 'saveahx':
        write_stream(main.file_id[:7]+'.ahx', main.st, overwrite=False)
        main.recalc = False

    elif uin == 'print station':
        print_station(main)
        main.recalc = False

    elif uin == 'pick segments' or uin == 'p':
        main = set_pick(main.spec, main, main.i, kwargs, set_all=False,
                        set_following=False)

    elif uin == 'replace pick' or uin == 'rp':
        main = replace_pick(main.spec, main, main.i, main.j,
                            kwargs, set_all=False,
                            set_following=False)
        main.update_maps = False

    elif uin == 'pick segments all' or uin == 'pa':
        main = set_pick(main.spec, main, main.i, kwargs, set_all=True,
                        set_following=False)
        main.update_maps = False

    elif uin == 'pick segments follow' or uin == 'pf':
        main = set_pick(main.spec, main, main.i, kwargs, set_all=False,
                        set_following=True)
        main.update_maps = False

    elif uin == 'print segments':
        if None in kwargs:
            sfile = "%s.pickle" % main.file_id
        else:
            sfile = "%s.%s" % (main.file_id, kwargs[0])
        main.segments.write(sfile, overwrite=True)
        main.recalc = False

    elif uin == 'printfpeaks' or uin == 'fpeaks' or uin == 'print fpeaks':
        main.st_work = main.st.copy()
        main.stsyn_work = main.st_syn[:]
        prep_fpeaks(main)
        main.recalc = False

    elif uin == 'load segments' or uin == 'load segment':
        main.segments = load_segments(kwargs)
        main.update_maps = True

    elif uin == 'unload segments':
        msg = "segment-file unloaded, no picks are shown anymore."
        print(msg)
        main.segments = None
        main.update_maps = False

    elif uin == 'reset':
        msg = '\033[93mReloading original files\033[0m\n'
        if None in kwargs:
            msg += 'Delete picked segments? (y/n) -->'
            ans = input(msg)
        else:
            ans = kwargs[0]

        if ans in ['y', 'Y']:
            if type(main.st) != list:
                main.st_work = main.st.copy()
                if main.st_syn is not None:
                    main.st_syn_work = main.st_syn[:]
            else:
                main.st_work = main.st[:]

            main.segments = None
            msg = '\033[93mSegments deleted, original files reloaded\033[0m\n'
            print(msg)
            main.save_stream = False
        elif ans in ['n', 'N']:
            if type(main.st) != list:
                main.st_work = main.st.copy()
                if main.st_syn is not None:
                    main.st_syn_work = main.st_syn[:]
            else:
                main.st_work = main.st[:]

            if main.segment_file:
                main.segments = load_segments([main.segment_file])

            msg = '\033[93mOriginal files reloaded\033[0m\n'
            print(msg)
            main.save_stream = False
        else:
            msg = '\033[93mReset aborted\033[0m\n'
            print(msg)

    elif uin == 'reload':
        main = reload(main, kwargs)

    elif uin == 'qcycle':
        qc, qwindow, modes = qcycle(main.modes_all, kwargs)
        if modes:
            main.fw[0] = modes[0].freq - 0.15
            main.fw[1] = modes[0].freq + 0.15
        main.update_maps = False

    elif uin == 'taper':
        tshape = set_taper()
        if tshape is None:
            main.recalc = False
        else:
            main.taper_shape = tshape
        main.update_maps = False

    elif uin == 'taper test':
        call_taper_test(main.tr, main.fw[0], main.fw[1], main.tw, main.modes,
                        kwargs)

    elif uin == 'tw test end':
        main.tw_old = main.tw[:]
        main.tw = init_time_window_test(main, kwargs, vpoint='end')

        if main.tw == main.tw_old:
            main.recalc = False
        else:
            main.tw_org = main.tw[:]

    elif uin == 'tw test start':
        main.tw_old = main.tw[:]
        main.tw = init_time_window_test(main, kwargs, vpoint='start')

        if main.tw == main.tw_old:
            main.recalc = False
        else:
            main.tw_org = main.tw[:]

    elif uin == 'font size' or uin == 'fs':
        main.fs = get_font_size(kwargs)

    # elif uin == 'nwin':
    #     if pick:
    #         if None in kwargs:
    #             nwin[-1] = True
    #         elif kwargs[0] == 'off':
    #             nwin[-1] = False
    #         elif kwargs[0] is False:
    #             nwin[-1] = False
    #         elif kwargs[0] in ['all', 'All']:
    #             nwin[-1] = True
    #             plot_noisewin(rfig, rax)
    #     elif nwin[-1] is True:
    #         if kwargs[0] == 'off':
    #             nwin[-1] = False
    #         elif kwargs[0] is False:
    #             nwin[-1] = False
    #     else:
    #         plot_noisewin(rfig, rax)

    elif uin == 'autopick' or uin == 'ap':

        main.segments = init_auto_pick(main.spec, main, kwargs)[0]

    elif uin == 'replace autopick' or uin == 'rap':

        main = replace_auto_pick(main.spec, main, main.j, kwargs)

    # elif uin == 'autopick all' or uin == 'apa':
    #     try:
    #         ifw, itw, min_snr, max_mf = init_autopick_loop(kwargs)
    #         main.segments = spectrum_loop(st_work=main.st_work,
    #                                       st_syn_work=main.st_syn_work,
    #                                       segments=main.segments, ifw=ifw,
    #                                       itw=itw,
    #                                       cmt=main.cmt,
    #                                       itaper_shape=main.taper_shape,
    #                                       iweighting=main.weighting,
    #                                       imin_snr=min_snr,
    #                                       imax_mf=max_mf)
    #     except Exception:
    #         msg += "\nError message:\n%s" % sys.exc_info()
    #         print(msg)

    elif uin == 'update segments' or uin == 'us':
        main.segments = update_segments(kwargs, None, main)
    elif uin == 'peak hw':
        init_twin_halfwidth(main.spec, main.modes, main.fw, kwargs)

    # elif uin == 'lat plot':
    #     init_plot_speclat(st_work, tr, [Htstart, Htend], [wstart, wend],
    #                       segments, modedata, kwargs)
    elif uin == 'zoom':
        main.zoom_sl = set_zoom_sl(main.zoom_sl, kwargs)

    elif uin == 'cmap' or uin == 'colormap':
        if None not in kwargs:
            if hasattr(cm, kwargs[0]):
                main.cmap = kwargs[0]
            else:
                msg = 'Unkown colormap\nSee: '
                msg += '(https://matplotlib.org/users/colormaps.html)'
                print(msg)

    else:
        main.recalc = False

    return pick, main


def set_zoom_sl(zoom_sl, kwargs):

    if None in kwargs:
        return zoom_sl
    else:
        return float(kwargs[0])


def replace_auto_pick(spec, main, j, kwargs):

    pick_old = main.segments.select(station=spec.stats.station.code)[j]
    main.segments = init_auto_pick(spec, main, kwargs)[0]
    main.segments.remove(pick_old)

    return main


def replace_pick(spec, main, i, j, kwargs, set_all, set_following):
    """
    Return the segment with the last pick replaced
    """

    segments = main.segments
    pick_old = segments.select(station=spec.stats.station.code)[j]
    main = set_pick(spec, main, i, kwargs, set_all, set_following)

    main.segments.remove(pick_old)
    return main


def reload(main, kwargs):
    return main


def init_twin_halfwidth(spec, modes_all, fw, kwargs):

    verbose = False
    interp = False
    max_peak = 'data'
    modes = []
    try:
        if None in kwargs:
            msg = 'Which mode? (name, e.g. 0S15) -->  '
            ans = input(msg)
            modes = get_mode(ans, modes_all)

            msg = 'Add extra plot, containing peak and range? (y,n) -->  '
            ans = input(msg)
            if ans in ['y', 'Y', 'yes', 'Yes']:
                verbose = True

            msg = 'Create envelope based on peaks in Amp.-Spectrum?'
            msg += ' (y,n) -->  '
            ans = input(msg)
            if ans in ['y', 'Y', 'yes', 'Yes']:
                interp = True

                msg = 'maximum peak based on interpolated data? (y,n) -->  '
                ans = input(msg)
                if ans in ['y', 'Y', 'yes', 'Yes']:
                    max_peak = 'interpolated'
        else:
            for arg in kwargs:
                if 'verbose' in arg:
                    verbose = True
                elif 'interpolate' in arg:
                    interp = True
                else:
                    modes = arg
        mode = get_mode(modes, modes_all)

        # try:
        hw, max_amp = calc_twin_halfwidth(spec, fw, mode,
                                          max_peak=max_peak,
                                          interpolate=interp,
                                          verbose=verbose)
        print("\n\tHalf-Width:  %.3e\n\tMax Amp:  %.3e" % (hw, max_amp))
        print("\n\tAHR:  %.3e" % (max_amp/hw))
    except Exception:
        msg = "Couldn't perform timewindow test. \n"
        msg += "Error message: %s" % sys.exc_info()[1]
        print(msg)
    return


def inspect_trace(tr, st_org, st, nstart, nend, t, i):
    tr_work = tr.copy()
    station = tr.stats.station
    title = "Pick Timewindow to be cut out"
    xlabel = 'time (h)'
    ylabel = ''

    x_axis = range(nstart, nend+1)
    xticks = np.linspace(nstart, nend, 11).astype(int)
    xtick_label = np.array(t[xticks] / 3600).astype(int)
    y_axis = tr.data[nstart:nend+1]
    picks = pick_segment_window(x_axis, y_axis, xticks=xticks,
                                xtick_label=xtick_label,
                                xlabel=xlabel, ylabel=ylabel,
                                title=title, pairwisepick=False)

    if len(picks) != 0:
        try:
            for win in pairwise(picks):
                start_index = win[0]
                end_index = win[1]
                tr_work.data = mask_data(tr_work.data, start_index, end_index,
                                         shape='linear')
            st.select(station=station)[0].data = tr_work.data
            st_org.select(station=station)[0].data = tr_work.data
            tr = tr_work
            msg = 'Data removed'
            save_stream = True
        except ValueError:
            save_stream = False
            st[i] = tr
        except Exception:
            save_stream = False
            msg = 'An error occured'
            st[i] = tr
    else:
        save_stream = False
        msg = 'No window selected'
        st[i] = tr

    print(msg)
    return tr, st_org, st, save_stream


def goto(st, i, kwargs):
    """
    1. jumps to station with input index
    2. jumps to first occurence of station with matching input name
    """
    if None in kwargs:
        station = input('Enter name/index of station  > ')
    else:
        station = kwargs[0]

    if station.isdigit():
        return int(station)
    else:
        for _i, tr in enumerate(st):
            if tr.stats.station == station:
                i = _i
        return i


def plot_noisewin(rfig, rax):
    return


def get_font_size(kwargs):
    if None in kwargs:
        fs = input('Enter new font size ->  ')
    else:
        fs = kwargs[0]
    return fs


def set_taper():
    msg = 'Set taper shape (Hanning, Hamming, Blackman, Kaiser, Bartlett)'
    msg += ' -->  '
    ans = input(msg)
    ans = ans.lower()
    if ans in ['hanning', 'hamming', 'blackman', 'kaiser', 'bartlett']:
        msg = 'Setting taper shape to %s' % ans[0].upper() + ans[1:]
        print(msg)
    else:
        msg = 'Wrong taper! (Hanning, Hamming, Blackman, Kaiser or Bartlett)\n'
        print(msg)
        ans = None

    return ans


def qcycle(modedata, kwargs=[None]):
    verbose = False
    if None in kwargs:
        ans = input('Which mode? (name, e.g. 0S15) -->  ')
        mode = get_mode(ans, modedata)
        verbose = True
    else:
        for kw in kwargs:
            if kw == 'verbose':
                verbose = True

        mode = get_mode(kwargs[0], modedata)

    if mode is None:
        return

    if verbose:
        msg = '\nQ-cycle length is %.1f h' % mode[0].qcycle
        print(msg)

    qwindow = (1, mode[0].qcycle)
    return mode[0].qcycle, qwindow, mode


def set_pick(spec, main, i, kwargs=None, set_all=False, set_following=False):

    st = main.st
    segments = main.segments
    fw = main.fw
    weighting = main.weighting
    modes = main.modes

    # Index range for freq axis:
    f_ind = np.arange(spec.flabel(fw[0]), spec.flabel(fw[1])+1)

    amp = abs(spec.data.fft.Data)[f_ind]
    if len(spec.syn) != 0:
        ampsyn = abs(list(spec.syn[0].fft.values())[0])[f_ind]
    else:
        ampsyn = None

    freq = spec.stats.freq[f_ind]
    xlabel = 'frequency (mHz)'
    ylabel = 'Amplitude'

    if None not in kwargs:
        seg_picks = [0, 0]
        try:
            seg_picks[0] = float(kwargs[0])
            seg_picks[1] = float(kwargs[1])
        except Exception:
            msg = 'Usage: p - fstart fend'
            print(msg)
            while True:
                try:
                    seg_picks[0] = input('Input fstart in mHz > ')
                    seg_picks[1] = input('Input fend   in mHz > ')
                    break
                except SyntaxError:
                    continue

    else:
        seg_picks = pick_segment_window(freq, amp, ampsyn,
                                        title='Mark 2 picks',
                                        xlabel=xlabel, ylabel=ylabel,
                                        modes=modes)
    seg_picks.sort()
    if len(seg_picks) != 0:
        segments = append_segment(spec, seg_picks, segments, weighting,
                                  cmt=main.cmt_id)
        if set_following:
            msg = 'Setting windows for all following stations'
            print(msg)
            tw = spec.stats.tw / 3600.
            for trace in st[i:]:
                print("%5s" % trace.stats.station, end='\r')
                spec = Spectrum(trace, tw[0], tw[1])
                segments = append_segment(spec, seg_picks, segments, weighting,
                                          cmt=main.cmt_id)
        if set_all:
            msg = 'Setting windows for all stations'
            print(msg)
            tw = spec.stats.tw / 3600.
            for trace in st:
                print("%5s" % trace.stats.station, end='\r')
                spec = Spectrum(trace, tw[0], tw[1])
                segments = append_segment(spec, seg_picks, segments, weighting,
                                          cmt=main.cmt_id)

    main.segments = segments
    return main


def init_auto_pick(spec, main, kwargs):

    modes_all, tw, tw_org = main.modes_all, main.tw, main.tw_org
    taper_shape, weighting = main.taper_shape, main.weighting
    seg = main.segments

    verbose = False
    min_snr = None
    max_mf = None
    modes = []
    if None in kwargs:
        ans = input('Which mode? (name, e.g. 0S15) -->  ')
        modes = get_mode([ans], modes_all)
        min_snr = input('Which minimum SNR? -->  ')
        if min_snr in ['none', 'None']:
            min_snr = None
        else:
            min_snr = float(min_snr)
        max_mf = input('Which maximum Misfit? -->  ')

        if max_mf in ['none', 'None']:
            max_mf = None
        else:
            max_mf = float(max_mf)

        ans = input('Add extra plot, containing peak and range? (y,n)')
        if ans in ['y', 'Y', 'yes', 'Yes']:
            verbose = True
    else:
        for arg in kwargs:
            if 'snr' in arg:
                min_snr = float(arg[-1])
            elif 'mf' in arg:
                max_mf = float(arg[-1])
            elif 'verbose' in arg:
                verbose = True
            else:
                modes.append(arg)
        modes = get_mode(modes, modes_all)

    segments = auto_set_pick(spec, tw, tw_org, seg, taper_shape, weighting,
                             modes, verbose)
    return segments, min_snr, max_mf, verbose, modes


def init_autopick_loop(kwargs):
    try:
        min_snr = None
        max_mf = None
        tw = [float(kwargs[0]), float(kwargs[1])]
        fw = [float(kwargs[2]), float(kwargs[3])]

        if len(kwargs) == 5:
            if kwargs[4] in ['none', 'None']:
                min_snr = None
            else:
                min_snr = float(kwargs[4])
        elif len(kwargs) == 6:
            if kwargs[5] in ['none', 'None']:
                max_mf = None
            else:
                max_mf = float(kwargs[5])

        return fw, tw, min_snr, max_mf

    except Exception:
        print('Usage: apa - tw1 tw2 fw1 fw2 (min_snr max_mf)')
        return


def auto_picker(spectrum, modes):

    fw1 = []
    fw2 = []
    flabels = []
    modelabels = []
    Fxx = spectrum.data.fft.Data
    f = spectrum.stats.freq
    if modes is not None:
        for mode in modes:
            _f = mode.freq
            # Define range
            delta = 0.02
            frange = [_f - delta, _f + delta]
            modelabel = get_freq_index(_f, spectrum.stats.delomeg)
            startlabel = get_freq_index(frange[0], spectrum.stats.delomeg)
            endlabel = get_freq_index(frange[1], spectrum.stats.delomeg)

            modelabels.append(modelabel)
            flabels.append([startlabel, endlabel])
    else:
        flabels = [[startlabel, endlabel]]
        modelabels = [int(np.round((endlabel-startlabel)/2. + startlabel))]

    for modelabel, flabel in zip(modelabels, flabels):
        startlabel = flabel[0]
        endlabel = flabel[1]
        # Find peaks
        spec = abs(Fxx[startlabel:endlabel+1])
        close_label = modelabel-startlabel
        peak_ind = find_peaks(spec, closest=close_label)[0]
        if not peak_ind:
            continue

        # Find neighbouring local minima
        minima = neighbouring_minima(abs(Fxx), startlabel+peak_ind[0])
        if None not in minima:
            fw1.append(f[startlabel+peak_ind[0] + minima[0]])
            fw2.append(f[startlabel+peak_ind[0] + minima[1]])
    if len(fw1) != 0:
        fw1 = np.array(fw1).min()
        fw2 = np.array(fw2).max()
        return fw1, fw2
    else:
        return None


def auto_set_pick(spec, tw, tw_org, seg, taper_shape, weighting, modes,
                  verbose=False, single_seg=False):

    # Write frequencies of modes
    freq = []
    for mode in modes:
        freq.append(mode.freq)

    f = spec.stats.freq
    Fxx = spec.data.fft.Data
    fw1 = []
    fw2 = []
    if verbose:
        fig, ax = plt.subplots()
    for _f in freq:
        # Define range
        delta = 0.02
        frange = [_f - delta, _f + delta]

        modelabel = spec.flabel(_f) + 1
        startlabel = spec.flabel(frange[0])
        endlabel = spec.flabel(frange[1])+1

        # Find peaks
        amp = abs(Fxx[startlabel:endlabel+1])
        close_label = modelabel-startlabel
        peak_ind = find_peaks(amp, closest=close_label)[0]

        if not peak_ind:
            continue
        # Find neighbouring local minima
        leftfound = False
        rightfound = False
        minima = [
                  abs(Fxx[startlabel+peak_ind[0]]),
                  abs(Fxx[startlabel+peak_ind[0]])
                  ]
        _i = 0
        while True:
            ri = _i
            x = abs(Fxx[startlabel+peak_ind[0] + _i])
            xx = abs(Fxx[startlabel+peak_ind[0] + ri + 1])
            if not rightfound:
                if x <= xx:
                    minima[1] = int(ri)
                    rightfound = True

            li = - _i
            y = abs(Fxx[startlabel+peak_ind[0] + li])
            yy = abs(Fxx[startlabel+peak_ind[0] + li - 1])
            if not leftfound:
                if y <= yy:
                    minima[0] = int(li)
                    leftfound = True

            if leftfound and rightfound:
                break
            _i += 1

        fw1.append(f[startlabel+peak_ind[0] + minima[0]])
        fw2.append(f[startlabel+peak_ind[0] + minima[1]])

    # Plot peaks
        if verbose:
            ax.plot(f[startlabel:endlabel+1], abs(Fxx[startlabel:endlabel+1]))
            ax.axvline(f[startlabel+peak_ind[0]])
            ax.axvspan(f[startlabel+peak_ind[0] + minima[0]],
                       f[startlabel+peak_ind[0] + minima[1]], alpha=0.2,
                       color='orange')

            ax.plot(f[startlabel+peak_ind[0]],
                    abs(Fxx[startlabel+peak_ind[0]]),
                    marker='o')
            ax.set_ylabel('Amplitude')
            ax.set_xlabel('frequency (mHz)')

            ax = format_exponent(ax)

            trans = ax.get_xaxis_transform()

            for value in mode.values():
                fname = value.name
                freq = value.freq
                ax.text(freq, 1.03, r"$%s$" % fname, transform=trans,
                        rotation=45, ha="center", va='center')
                ax.axvline(x=freq, linestyle=':', linewidth=1, color='grey')

    if fw1 and fw2:
        fw1 = np.array(fw1).min()
        fw2 = np.array(fw2).max()
        segments = append_segment(spec, [fw1, fw2], seg, weighting)
    else:
        segments = append_segment(spec, [frange[0], frange[1]], seg,
                                  weighting)

    if single_seg:
        segments = segments.select(station=spec.stats.station.code)

    return segments


def calc_twin_halfwidth(spec, fw, mode, max_peak='data', interpolate=True,
                        verbose=False):
    """
    Calculates the half-width around a given frequency peak in the amplitude
    spectrum. For this it checks for local maxima in the range from wstart
    to wend and interpolates the points using

    scipy.interp1d(kind='quadratic')

    The resulting interpolated function is the used to calculated the half
    width of the frequency peak of the input mode.
    """

    startlabel = spec.flabel(fw[0])
    endlabel = spec.flabel(fw[1])
    amp = abs(spec.data.fft.Data[startlabel:endlabel+1])
    a = get_local_extrema(amp[startlabel:endlabel+1], 'max')
    a += startlabel

    if not interpolate or len(a) < 3 or SCIPY is False:
        env = amp[startlabel:endlabel+1].copy()
        xf1 = np.arange(startlabel, endlabel+1)
        title = 'Full Width Half Magnitude for data'
    else:
        f1 = interp1d(a, amp[a], kind='quadratic')
        xf1 = np.arange(a.min(), a.max()+1)
        env = f1(xf1)
        title = 'Full Width Half Magnitude for interpolated data'
    # corresponding freqs:
    # f[xf1+startlabel]

    Amax = env.argmax()
    Ahw = np.where(env > env[Amax]/2.)[0]
    freq = []
    for name, value in mode.items():
        freq.append(value.freq)

    modelabel = spec.flabel(freq[0]) + 1 - startlabel

    before, after = takeClosest(Ahw, modelabel, N=2)

    # If modelabel is in AHw, after will have the same value
    if after == modelabel:
        before = after
    cons_Ahw = group_consecutive_numbers(Ahw)
    hw_points = [0, 0]
    for c in cons_Ahw:
        if before in c or after in c:
            if env[c].max() > hw_points[0]:
                hw_points = [env[c].max(), np.array(c)]

    if max_peak == 'data':
        ymax = amp[xf1[0] + hw_points[1]].max()
        xmax = amp[xf1[0] + hw_points[1]].argmax()+xf1[0]+hw_points[1][0]
    elif max_peak == 'interpolated':
        xmax = xf1[0] + hw_points[1][0] + env[hw_points[1]].argmax()
        ymax = hw_points[0]

    if verbose:
        f = spec.stats.freq
        fig, ax = plt.subplots()
        ax.plot(f[startlabel:endlabel+1], amp[startlabel:endlabel+1])
        ax.plot(f[xf1], env, 'grey')
        # Ahw_sl = Ahw + a.min()
        # ax.plot(f[Ahw_sl], env[Ahw], color='red')

        ax.plot(f[xf1[0] + hw_points[1]], env[hw_points[1]], 'r--',
                linewidth=3)

        xmax = xf1[0] + hw_points[1][0] + env[hw_points[1]].argmax()

        ax.plot(f[xmax], ymax, 'ro')

        ax.hlines(env[Amax]/2., f[xf1[0] + hw_points[1][0]],
                  f[xf1[0] + hw_points[1][-1]])

        ax.plot(f[xf1[0] + hw_points[1][0]], env[Amax]/2., 'b|', linewidth=5)
        ax.plot(f[xf1[0] + hw_points[1][-1]], env[Amax]/2., 'b|', linewidth=5)
        ax.set_ylabel('Amplitude')
        ax.set_xlabel('frequency (mHz)')
        ax.set_title(title)

    hw = abs(f[xf1[0] + hw_points[1][0]] - f[xf1[0] + hw_points[1][-1]])

    return hw, ymax


def read_singlets(trace):
    sname = "*" + trace.stats.station + "*singlets"
    singfile = str(glob.glob(sname)[0])
    singlets = np.genfromtxt(singfile)
    if len(singlets) == 2:
        singlets = np.array([singlets])
    # singlets = singlets.transpose()

    return singlets, singfile


def read_omega(omegafile='omega.dat'):
    omega = np.genfromtxt(omegafile)
    if len(omega) == 2:
        omega = np.array([omega])
    return omega


def remove_station(main, save_tmp=False, all=False):

    if all is True:
        st_tmp = main.st
    else:
        st_tmp = main.st.select(station=main.tr.stats.station)

    for tr in st_tmp:
        main.st.remove(tr)
        stw = main.st_work.select(station=tr.stats.station)
        for trw in stw:
            main.st_work.remove(trw)

        if main.segments:
            try:
                ps = main.segments.select(station=main.tr.stats.station)
                for p in ps:
                    main.segments.remove(p)
            except ValueError:
                pass

    if all is False:
        msg = '\n\033[93mStation: %s has ' % main.tr.stats.station
        msg += 'been removed\033[0m'
    else:
        msg = '\n\033[93mAll stations have been removed\033[0m'

    print(msg)
    if save_tmp:
        name = '%s.tmp.ahx' % main.file_id
        save_data(main, name=name, overwrite=True)
    if main.st_syn:
        for sts in main.st_syn:
            st_tmp = sts.select(station=main.tr.stats.station)
            for tr in st_tmp:
                sts.remove(tr)

        for sts in main.st_syn_work:
            st_tmp = sts.select(station=main.tr.stats.station)
            for tr in st_tmp:
                sts.remove(tr)
        if save_tmp:
            for snum, sts in enumerate(main.st_syn):
                name = '%s.%s.tmp.ahx.syn' % (main.file_id, str(snum))
                save_data(main, name)
    return


def remove_pick(spec, segments, kwargs, verbose=True):
    picks = segments.select(station=spec.stats.station.code)
    if len(kwargs) > 1:
        msg = ''
        for p_nums in kwargs[1:]:
            p = picks[int(p_nums)]
            segments.remove(p)
            msg += '\n\033[93mPick: %s has been removed\033[0m' % p

    elif len(picks) > 1:
        msg = "More than 1 pick to delete for station %s"
        msg = msg % spec.stats.station.code
        for i, p in enumerate(picks):
            msg += "\n%i: %s" % (i, p)

        msg += '\n\nEnter pick number (sep. by comma)\n'
        msg += 'Waiting for user input\n -->  '
        uinput = input(msg)

        msg = ''
        if uinput.isdigit():
            p = picks[int(uinput)]
            segments.remove(p)
            msg += '\n\033[93mPick: %s has been removed\033[0m' % p
        elif len(uinput) > 1:
            for p_nums in uinput.split(','):
                p = picks[int(p_nums)]
                segments.remove(p)
                msg += '\n\033[93mPick: %s has been removed\033[0m' % p
        else:
            return

    else:
        segments.remove(picks[0])
        msg = '\n\033[93mPick: %s has been removed\033[0m' % picks[0]

    if verbose:
        print(msg)
    return


def save_data(main, name=None, overwrite=False):
    # if not name:
    #     name = '%s.ahx' % file_id
    # write_stream(name, st, overwrite)
    if main.segments:
        segfile = "%s.pickle" % main.file_id
        main.segments.write(segfile, overwrite)
    return


def setup_plot_spectrum(gs, n_plots, minispec):
    fig, ax = plt.subplots()
    ax = range(n_plots)
    ax = AttribDict()
    if minispec:
        # Timeseries
        ax['seis'] = plt.subplot(gs[0:2, 0])
        # Amplitudes
        ax['amp'] = plt.subplot(gs[3:5, 0])
        seg_ax = [ax.amp]
    else:
        # # Real part
        ax['real'] = plt.subplot(gs[0:2, 0:2])

        # # Imaginary part
        ax['imag'] = plt.subplot(gs[2:4, 0:2])

        # Phase
        ax['phase'] = plt.subplot(gs[4:6, 0:2])

        # Amplitude
        ax['amp'] = plt.subplot(gs[6:12, 0:2])

        # Timeseries
        ax['seis'] = plt.subplot(gs[10:12, 2:4])

        # Magnitudes
        ax['mag'] = plt.subplot(gs[8:10, 2:4])

        # These axis is set in maps
        # ax[4] = init_gcpmap(gs[4:8, 2:4], slon, elon)
        #
        # if inv:
        #     ax[5] = init_map(gs[0:4, 2:4])
        # Latitude plot of st_work
        ax['speclat'] = plt.subplot(gs[0:12, 4:6])
        seg_ax = [ax.real, ax.imag, ax.phase, ax.amp]
    return fig, ax, seg_ax


def get_iter_colormap(input_list, cmap):
    if len(input_list) > 3:
        colormap = iter(getattr(cm, cmap)(
                        np.linspace(0, 1, len(input_list)))
                        )
    else:
        colormap = iter(['blue', 'red', 'green'])
    return colormap


def init_map(gs):
    """
    param grd_size: defining the grid size
    type  grd_size: 2x2 tupel

    param grd_pos: position on grid
    type  grd_pos: 2x2 tupel
    """
    ocean_color = '#EBEBEB'
    land_color = '#FBFBF2'
    ax = plt.subplot(gs, projection=ccrs.Mollweide())
    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()
    # ax.stock_img()
    ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    ax.coastlines()
    ax.gridlines()
    return ax


def init_gcpmap(gs, lon1, lon2):
    ocean_color = '#EBEBEB'
    land_color = '#FBFBF2'

    difflon = abs(lon1 - lon2)
    if difflon > 180:
        lon0 = int((lon1+lon2)/2. + 180)
    else:
        lon0 = int((lon1+lon2)/2.)
    proj = ccrs.Mollweide(central_longitude=lon0)
    ax = plt.subplot(gs, projection=proj)
    ax.set_global()
    # ax.stock_img()
    ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    ax.coastlines()
    ax.gridlines()
    return ax


def _plot_modes_in_ax(ax, modes, f_axis, color='grey'):
    trans = ax.get_xaxis_transform()
    for mode in modes:
        if mode.freq > f_axis[0] and mode.freq < f_axis[1]:
            ax.text(mode.freq, 1.03, r"$%s$" % mode.name,
                    transform=trans, rotation=45, ha="center",
                    va='center')
            ax.axvline(x=mode.freq, linestyle=':', linewidth=1, color=color)
    return


def update_segments(kwargs, spec, main):
    """
    Recalculate the weights for a given segment file
    """
    if None in kwargs:
        msg = "\033[93mGive kwarg: tw, fw, or weighting\033[0m"
        print(msg)
        return main.segments
    else:
        for arg in kwargs:
            if 'tw' in arg:
                recalc = 'tw'
            elif 'fw' in arg:
                recalc = 'fw'
            elif 'weighting' in arg:
                recalc = 'weighting'
        msg = "\033[92mUpdating segment-file\033[0m"
        print(msg)

    if spec:
        seg = main.segments.select(station=spec.stats.station.code)
    else:
        seg = main.segments

    for i, picks in enumerate(seg):
        try:
            tr = main.st.select(station=picks.stats.station)[0]
        except IndexError:
            seg.remove(picks)
            continue

        if recalc != 'tw':
            main.tw = [picks.tw1, picks.tw2]
        # This might be unnecessairy!
        if recalc != 'fw':
            old_fw = main.fw
            main.fw = [picks.fw1, picks.fw2]

        # Check if spectrum object is the same as station -> no FFT calc
        if spec is None:
            spec = Spectrum(tr, main.tw[0], main.tw[1], taper=main.taper_shape)
        elif tr.stats.station != spec.stats.station.code:
            spec = Spectrum(tr, main.tw[0], main.tw[1], taper=main.taper_shape)
        elif recalc == 'tw':
            spec = Spectrum(tr, main.tw[0], main.tw[1], taper=main.taper_shape)

        main.segments.remove(picks)
        main.segments = append_segment(spec, seg_picks=main.fw,
                                       segments=main.segments,
                                       weighting=main.weighting,
                                       cmt=main.cmt)
        main.segments.add_channel(tr.stats.channel)

        main.fw = old_fw

    return main.segments


def init_plot_speclat(st, tr, tw, fw, segments, modedata, kwargs):
    zoom = 1
    color = 'black'
    highlight = tr.stats.station
    ax = None
    tr_ref = None
    #
    if None not in kwargs:
        for kw in kwargs:
            if kw.isdigit():
                zoom = kw

    plot_speclat(st, tw, fw, zoom, segments, tr_ref, highlight, modedata,
                 ax, color)

    return


def get_freq_index(freq, delomeg):
    return int(np.round(freq * 2. * np.pi/(1000. * delomeg))) + 1
