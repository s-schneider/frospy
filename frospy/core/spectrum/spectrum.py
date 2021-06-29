# -*- coding: utf-8 -*-
"""
Module for handling nmPy segment objects.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import matplotlib.pyplot as plt
from frospy.plot.nmplt import (get_iter_colormap, get_iter_linestyle,
                               format_exponent)

from frospy.core.segment import Segment
from frospy.util.base import (inv4trace, nextpow2, freq_index, time_index,
                              get_local_extrema, group_consecutive_numbers,
                              takeClosest, misfit)
from frospy.util.base import taper as taper_trace
from frospy import data as frospydata

import numpy as np
from obspy.core import AttribDict
from obspy.core.inventory import Station
from obspy.core.trace import Stats as TraceStats
# from copy import deepcopy

try:
    from scipy.interpolate import interp1d
    SCIPY = True
except ImportError:
    SCIPY = False


class Stats(AttribDict):
    """
    Stats object for Spectrum class
    """
    defaults = {
        'tw': np.array([]),
        'taper': 'hanning',
        'station': None,
        'record': None,
        'delomeg': 0,
        'freq': None,
        'origin_Tdiff': 0,
        'times': None
                }

    _refresh_keys = {'tw', 'taper', 'station', 'record', 'delomeg', 'freq',
                     'origin_Tdiff', 'times'}

    def __init__(self, header={}):
        """
        """
        # super() is used to initialize AttribDict within Stats(),
        # which sets the defaults as initial dictionary
        # equivalent to AttribDict.__init__(self, header)
        super(Stats, self).__init__(header)

    def __setitem__(self, key, value):
        """
        """
        if key in self._refresh_keys:
            if key in ['taper', 'hanning']:
                value = str(value)
            if key in ['tw']:
                if type(value) is list:
                    value = value
            elif key in ['delomeg', 'origin_Tdiff']:
                value = float(value)
            elif key in ['freq', 'times']:
                if isinstance(value, np.ndarray):
                    value = value
            elif key == 'station':
                if isinstance(value, Station):
                    value = value
            elif key == 'record':
                if isinstance(value, TraceStats):
                    value = value
            # equivalent to AttribDict.__setitem__(self, key, value)
            super(Stats, self).__setitem__(key, value)

    __setattr__ = __setitem__

    def __str__(self, extended=False):
        """
        Return better readable string representation of Stats object
        """
        _pretty_str = '%s | %s | %s-%s hrs' % (self.station.code,
                                               self.record.channel,
                                               self.tw[0]/3600.,
                                               self.tw[1]/3600.)
        return _pretty_str

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))


class Data(object):
    """
    Data object for Spectrum class
    """
    # calc FFT here and set additional self.stats values,
    # e.g. self.stats.delomeg
    def __init__(self, trace, stats, label):
        t1 = stats.tw[0]
        t2 = stats.tw[1]

        power = nextpow2(len(trace.data))
        if power < 16:
            power = 16
        maxdata = np.power(2, power)

        mid = t1 + 0.5*(t2 - t1)
        shift_fac = np.exp(1j * stats.freq * mid * 2. * np.pi/1000.)

        tdiff = stats.origin_Tdiff
        nstart = int(round((t1-tdiff)/trace.stats.delta)+1)
        nend = int(round((t2-tdiff)/trace.stats.delta))
        tracefill = taper_trace(nstart, nend, trace, stats.taper)

        X = np.fft.fft(tracefill, maxdata)
        FT = X * shift_fac

        timeseries = AttribDict()
        timeseries[label] = tracefill
        fft = AttribDict()
        fft[label] = FT

        self.timeseries = timeseries
        self.fft = fft

    def __str__(self, extended=False):
        for key in self.timeseries.keys():
            out = "Spectrum for: %s\n" % key
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))


class Spectrum(object):
    """
    Class that handles the spectrum produced by
    frospy.preprocessing.spectrum methods

    """

    def __init__(self, trace, tw_start, tw_end, syn_stream=None, label=None,
                 syn_label=None, taper='hanning'):
        """
        trace: trace object, containing data of one station
        tw_start, tw_end: starttime and endtime of timewindow in hours.
        syn_stream: stream object containing the corresponding synthetic
                    traces to 'trace'
        """
        # The order is important, first set Stats(header), then set Data
        header = create_header(trace, tw_start, tw_end, taper)
        self.stats = Stats(header)
        self.segments = Segment()

        # Add data
        if label is None:
            label = 'Data'
        self.data = Data(trace, self.stats, label)
        # Add synthetics
        if syn_stream is not None:

            if syn_label is not None and type(syn_label) == str:
                syn_label = [syn_label]

            if syn_label is not None and len(syn_stream) != len(syn_label):
                msg = 'Not enough syn_label provided!'
                print(msg)
                syn_label = None

            if syn_label is None:
                num = np.arange(1, len(syn_stream)+1).astype(str)
                name = num.copy()
                name[:] = 'synthetics '
                syn_label = np.core.defchararray.add(name, num)

            self.syn = []
            for i, syn_tr in enumerate(syn_stream):
                self.syn.append(Data(syn_tr, self.stats, syn_label[i]))
        else:
            self.syn = None

    def __str__(self, extended=False):
        out = "%s" % self.stats.__str__()
        out += "\n%s" % self.data.__str__()
        if self.syn is not None:
            for syn in self.syn:
                out += "%s" % syn.__str__()
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def weight(self, fw1, fw2, rule='sum'):
        nstart = freq_index(fw1, self.stats.delomeg)
        nend = freq_index(fw2, self.stats.delomeg)

        F = list(self.data.fft.values())[0]
        amp = abs(F[nstart:nend+1])
        if rule == 'integrate':
            return np.trapz(amp)
        else:
            return sum(amp)

    def signal2noise(self, pick):
        """
        Signal to noise ration, calculated by the ratio of the maximum value in
        the given frequency window for data (defined in pick.fw1 and pick.fw2).
        The noise value is calulated using the windows in
        nmpy/data/noisewindows.dat .

        It takes the closes windows to the left and right of the given
        frequencywindow and picks the maximum amplitude as the noise value

        All values are computed as the amplitude spectrum:
            abs(FFT) = sqrt(Re(FFT)^2 + Im(FFT)^2)

        param
        type

        returns
        """
        path = frospydata.__path__[0] + "/AD/noisewindows_new.dat"
        nwins = np.genfromtxt(path)

        # Finding correct noisewindows
        noise = None
        overlap_err = False
        for i, void in enumerate(nwins[0:-1]):
            if nwins[i][1] < pick.fw1:
                if nwins[i+1][0] > pick.fw2:
                    if noise is None:
                        noise = []
                    noise.append(nwins[i])
                    noise.append(nwins[i+1])
                    break

        # If none found, try to include data in noise window as well
        if noise is None:
            overlap_err = True
            # Loop for lower bound
            for i, void in enumerate(nwins[0:-1]):
                if nwins[i][1] < pick.fw1 and nwins[i+1][0] > pick.fw1:
                    if noise is None:
                        noise = []
                    noise.append(nwins[i])
                    break
                elif nwins[i][0] <= pick.fw1 and nwins[i][1] >= pick.fw1:
                    if noise is None:
                        noise = []
                    if i != 0:
                        noise.append(nwins[i-1])
                    else:
                        noise.append(nwins[0])
                    break

            # Loop for upper bound
            for i, void in enumerate(nwins[0:-1]):
                if (nwins[i][1] < pick.fw2 and nwins[i+1][0] > pick.fw2) or \
                   (nwins[i][0] <= pick.fw2 and nwins[i][1] >= pick.fw2):
                    if noise is None:
                        noise = []
                    noise.append(nwins[i+1])
                    break
        if noise is None:
            return 0, None, False

        ampnoise = []
        for n in noise:
            startlabel = freq_index(n[0], self.stats.delomeg)
            endlabel = freq_index(n[1], self.stats.delomeg)
            amp = abs(list(self.data.fft.values())[0][startlabel:endlabel+1])
            ampnoise.append(amp.max())

        ampnoise = np.array(ampnoise).max()

        startlabel = freq_index(pick.fw1, self.stats.delomeg)
        endlabel = freq_index(pick.fw2, self.stats.delomeg)
        ampsignal = abs(
                list(self.data.fft.values())[0][startlabel:endlabel+1]).max()

        return ampsignal/ampnoise, noise, overlap_err

    def misfit(self, pick, label='all'):
        """
        Calculates misfit between data and synthetics give in
        self.data.fft and self.syn.fft

        output: list [rmf, cmf]
        """
        startlabel = freq_index(pick.fw1, self.stats.delomeg)
        endlabel = freq_index(pick.fw2, self.stats.delomeg)
        X = list(self.data.fft.values())[0]

        mf = []
        for s in self.syn:
            if label != 'all' and label != list(s.fft.keys())[0]:
                continue
            Xsyn = list(s.fft.values())[0]
            mf.append(misfit(X, Xsyn, startlabel, endlabel))
        return mf

    def signal2fwhm(self, pick, max_peak='data', interpolate=True,
                    peak_order=1.):
        """
        Full width at half maximum for a peak at x_center
        Calculates the ratio:

        A_max / (fwhm * N-of-local-maxima^2)


        checks for the range of values given in y[1:-2]

        * : y
        | : calculation area, excluding first and last sample

           |  *  |
           | * * |
           |*   *|
          *|     |*
        ___|_____|___________

        """
        seg_ind1 = freq_index(pick.fw1, self.stats.delomeg)
        seg_ind2 = freq_index(pick.fw2, self.stats.delomeg)

        y = list(self.data.fft.values())[0][seg_ind1-1:seg_ind2+2]
        x_center = y.argmax()

        amp = abs(y[1:-1])
        amp_full = abs(y)
        a = get_local_extrema(amp, 'max')
        no_of_peaks = float(len(a))

        if not interpolate or no_of_peaks <= 4 or SCIPY is False:
            env = amp.copy()
            env_shift = 0
        else:
            f1 = interp1d(a, amp[a], kind='quadratic')
            xf1 = np.arange(a.min(), a.max()+1)
            env = f1(xf1)
            env_shift = a.min()

        Ahw = np.where(env >= env.max()/2.)[0]
        cons_Ahw = group_consecutive_numbers(Ahw)
        modelabel = x_center + 1

        # If condition that checks if Ahw is half open
        before, after = takeClosest(Ahw, modelabel, N=2)

        # If modelabel is in AHw, after will have the same value
        if after == modelabel or before is None:
            before = after

        hw_points = [0, 0]
        for c in cons_Ahw:
            if before in c or after in c:
                # Check at both ends if intervals are open
                #   -> set hw to inf if so
                # Check left end:
                if c[0] == 0:
                    if amp_full[0] > env.max()/2.:
                        return 0.0
                elif env[c[0]-1] > env.max()/2.:
                    return 0.0
                # Check Right end
                if c[-1] == len(env)-1:
                    if amp_full[-1] > env.max()/2.:
                        return 0.0
                elif env[c[-1] + 1] > env.max()/2.:
                    return 0.0

                if env[c].max() > hw_points[0]:
                    hw_points = [env[c].max(), np.array(c)]

        hw = abs(hw_points[0] - hw_points[-1])
        if hw == 0:
            hw = 1
            ymax = 0
        else:
            if max_peak == 'interpolated':
                ymax = env[hw_points[1]].max()
            elif max_peak == 'data':
                ymax = amp[env_shift + hw_points[1]].max()

        return ymax/hw/np.power(no_of_peaks, peak_order)

    def plot(self, fw1, fw2, part='Amplitude', ax=None, width=0.825,
             cmap='rainbow', xlabel='f(mHz)', ylabel=None, dlabel=None,
             normalize=False, ticks=None, cmap_highlight=None,
             color='k',
             **plotargs):

        if ax is None:
            fig, ax = plt.subplots()

        startlabel = freq_index(fw1, self.stats.delomeg)
        endlabel = freq_index(fw2, self.stats.delomeg)

        # Set freq axis
        f = self.stats.freq

        # Plot data
        Fxx = list(self.data.fft.values())[0]
        if dlabel is None:
            dlabel = list(self.data.fft.keys())[0]

        y = get_part(Fxx[startlabel:endlabel+1], part)
        if normalize is True:
            ynorm = max(y)
        else:
            ynorm = 1

        y = y / ynorm
        ax.plot(f[startlabel:endlabel+1], y, linestyle='solid', color=color,
                linewidth=width, label=dlabel, **plotargs)

        # Plot synthetics
        if self.syn is not None:
            colormap = get_iter_colormap(self.syn, cmap)
            linestyles = get_iter_linestyle(exclude='-')
            for _i, trsyn in enumerate(self.syn):
                fs = list(trsyn.fft.values())[0]
                if _i == 0 and cmap_highlight is not None:
                    c = cmap_highlight
                else:
                    c = next(colormap)
                y = get_part(fs[startlabel:endlabel+1], part)
                y = y / ynorm
                label = list(trsyn.fft.keys())[0]
                ax.plot(f[startlabel:endlabel+1], y,
                        linestyle=next(linestyles), linewidth=width, c=c,
                        label=label, **plotargs)

        # Labeling
        if part in ['Re', 'Im', 'Phase']:
            if ylabel:
                ax.set_ylabel(ylabel)
            else:
                ax.set_ylabel(part)
        else:
            ax.set_ylabel(ylabel)
        if xlabel is None:
            ax.get_xaxis().set_ticklabels([])
        else:
            ax.set_xlabel(xlabel)

        if normalize is True:
            ax.yaxis.set_ticks([0])
        if ticks:
            ax.yaxis.set_ticks(ticks)
        ax.set_xlim(f[startlabel], f[endlabel + 1])
        ax = format_exponent(ax)
        return ax

    def get_spectrum(self, fw1, fw2, part='Amplitude',
                     normalize=False):

        startlabel = freq_index(fw1, self.stats.delomeg)
        endlabel = freq_index(fw2, self.stats.delomeg)

        # Set freq axis
        f = self.stats.freq

        # Plot data
        Fxx = list(self.data.fft.values())[0]
        y = get_part(Fxx[startlabel:endlabel + 1], part)
        if normalize is True:
            ynorm = max(y)
        else:
            ynorm = 1

        y = y / ynorm
        return f[startlabel:endlabel + 1], y

    def flabel(self, f):
        return freq_index(f, self.stats.delomeg)

    def tlabel(self, t):
        return time_index((t-self.stats.origin_Tdiff), self.stats.record.delta)


def create_header(trace, tw_start, tw_end, taper):
    """
    tw_start, tw_end in hours
    tw in seconds
    """
    tw_start = tw_start * 3600.  # from hrs to seconds
    tw_end = tw_end * 3600.  # from hrs to seconds
    tw = time_sanity_check(trace, tw_start, tw_end)

    header = {}
    header.setdefault('tw', np.array(tw))
    header.setdefault('taper', taper)
    header.setdefault('station', inv4trace(trace)[0][0])
    header.setdefault('record', trace.stats)

    power = nextpow2(len(trace.data))
    if power < 16:
        power = 16
    maxdata = np.power(2, power)  # len(trace.data)
    indvec = np.arange(maxdata)
    delomeg = 2. * np.pi / (maxdata * trace.stats.delta)
    f = 1000. * indvec * delomeg / (2. * np.pi)

    header.setdefault('freq', f)
    header.setdefault('delomeg', delomeg)
    header.setdefault('origin_Tdiff', origin_time_diff(trace))

    t = np.arange(trace.stats.npts) * trace.stats.delta
    t = t + header['origin_Tdiff']
    header.setdefault('times', t)
    return header


def time_sanity_check(trace, tw_start, tw_end, verbose=False):
    # org input: tr, Htstart, Htend, Htstart_org, Htend_org, verbose=True
    # Set errors to false
    STshort = False
    STlong = False
    ETlong = False

    tr = trace
    tstart = origin_time_diff(tr)

    # if Htstart_org is not None:
    #     Htstart, Htend = Htstart_org, Htend_org

    nstart = time_index((tw_start-tstart), tr.stats.delta)
    nend = time_index((tw_end-tstart), tr.stats.delta) + 1
    t = np.arange(tr.stats.npts) * tr.stats.delta
    t = t + tstart

    if nstart < 0:
        Tdiff = tw_end - tw_start
        tw_start = tstart
        tw_end = tw_start + Tdiff
        nstart = 0
        nend = int(round((tw_end-tstart)/tr.stats.delta))
        t = t + tstart
        s = tw_start / 3600.
        e = tw_end / 3600.
        STshort = True

    elif nstart > len(tr.data):
        Tdiff = tr.stats.endtime - tr.stats.starttime
        tw_start = tr.stats.starttime - tr.stats.ah.event.origin_time
        tw_end = tw_start + Tdiff
        nstart = 0
        nend = int(round((tw_end-tstart)/tr.stats.delta))
        s = tw_start / 3600.
        e = tw_end / 3600.
        STlong = True

    if nend > tr.stats.npts:
        nend = tr.stats.npts - 1
        tw_end = np.floor(nend * tr.stats.delta + tstart)
        Tdiff = tw_end - tw_start
        ETlong = True

    if verbose:
        s = tw_start / 3600.
        e = tw_end / 3600.

        msg = ''
        if STshort or STlong or ETlong:
            msg = 'Station: %s' % tr.stats.station
        if STshort:
            msg += '\n\033[93mStarttime earlier then record time'

        if STlong:
            msg += '\n\033[93mStarttime later then record time'

        if ETlong:
            msg += '\n\033[93mEndtime later then record time'
        if msg and verbose:
            msg += '\nTimewindow set to %.1f-%.1f h\033[0m\n'
            print(msg % (s, e))
    return tw_start, tw_end


def origin_time_diff(trace):
    """
    Calculates the difference of record and origin starttime
    """
    if hasattr(trace.stats, 'ah'):
        rsec = trace.stats.starttime.second
        rsec = rsec + trace.stats.starttime.microsecond * 1E-6
        esec = trace.stats.ah.event.origin_time.second
        esec = esec + trace.stats.ah.event.origin_time.microsecond * 1E-6
        rmn = trace.stats.starttime.minute
        emn = trace.stats.ah.event.origin_time.minute
        rhr = trace.stats.starttime.hour
        ehr = trace.stats.ah.event.origin_time.hour
        nde = trace.stats.starttime.day
        ndr = trace.stats.ah.event.origin_time.day

        tdiff = rsec-esec+60.*((rmn-emn) + 60.*((rhr-ehr) + 24.*(nde-ndr)))

        return tdiff

    else:
        msg = 'Data not in AH format, no event information found'
        raise IOError(msg)


def get_part(Fxx, part='Amplitude'):
    if part == 'Re':
        return Fxx.real
    elif part == 'Im':
        return Fxx.imag
    elif part == 'Phase':
        return np.angle(Fxx)
    elif part == 'Amplitude':
        return abs(Fxx)
