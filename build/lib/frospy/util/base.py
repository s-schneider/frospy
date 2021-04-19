from __future__ import absolute_import
from bisect import bisect_left
from numbers import Number
import numpy as np
from numpy import diff, sign
import math
import obspy
import re
import sys
import os
from obspy.clients.fdsn import Client
from obspy.core.inventory import Inventory, Network, Station
from obspy import Stream, Trace
from obspy.core.util.attribdict import AttribDict
from obspy.core.event.base import QuantityError
from itertools import zip_longest, groupby
from operator import itemgetter
from natsort import natsorted

from frospy import data as frospydata

try:
    import scipy as sp
    from scipy.interpolate import interp1d
    SCIPY = True
except ImportError:
    SCIPY = False
"""
Basic collection of fundamental functions for the nmPy lib
Author: S. Schneider 2016
"""


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def format_list(input, format='exp', output='str'):
    if type(input) is not list:
        input = [input]

    outlist = []
    for d in input:
        if output == 'str':
            if format == 'exp':
                outlist += ["%4.2e" % float(d)]
            elif format == 'float':
                outlist += ["%s" % float(d)]
        elif output == 'float':
            outlist += [float(d)]

    return outlist


def listdir(dir, attrib, keys=None, hiddenfolder=False):
    """
    lists all directorys in dir, that start with attrib
    If keys is set, only matching strings are returned
    """
    dir_tmp = []
    for f in os.listdir(dir):
        if not hiddenfolder and f.startswith('.'):
            continue
        if not os.path.isdir('%s/%s' % (dir, f)):
            continue
        if f.startswith(attrib):
            dir_tmp.append(f)

    dirs = []
    if type(keys) == list:
        for i, d in enumerate(dir_tmp):
            if float(d.split(attrib)[-1]) in keys:
                dirs.append(d)
    else:
        dirs = dir_tmp

    return dirs


def misfit(X, Xsyn, s, e):
    d2 = sum(abs(X[s:e+1])**2.)
    rmf = sum((abs(X[s:e+1]) - abs(Xsyn[s:e+1]))**2.) / d2
    cmf = sum(abs(X[s:e+1] - Xsyn[s:e+1])**2.) / d2

    return rmf, cmf


def signal2noise(Fxx, pick, delomeg, snr_min=None):
    """
    Signal to noise ration, calculated by the ratio of the maximum value in
    the given frequency window for data (defined in pick.fw1 and pick.fw2).
    The noise value is calulated using the windows in
    frospy/data/noisewindows_new.dat .

    It takes the closes windows to the left and right of the given frequency
    window and picks the maximum amplitude as the noise value

    All values are computed as the amplitude spectrum:
        abs(FFT) = sqrt(Re(FFT)^2 + Im(FFT)^2)

    param
    type

    returns
    """
    path = frospydata.__path__[0] + "/AD/noisewindows_new.dat"
    nwins = np.genfromtxt(path)

    # Finding correct noisewindows
    noise = []
    nind = []
    overlap_err = False
    for i, void in enumerate(nwins[0:-1]):
        if nwins[i][1] < pick.fw1:
            if nwins[i+1][0] > pick.fw2:
                noise.append(nwins[i])
                noise.append(nwins[i+1])
                break

    if not noise:
        overlap_err = True
        # Loop for lower bound
        for i, void in enumerate(nwins[0:-1]):
            if nwins[i][1] < pick.fw1 and nwins[i+1][0] > pick.fw1:
                noise.append(nwins[i])
                break
            elif nwins[i][0] <= pick.fw1 and nwins[i][1] >= pick.fw1:
                noise.append(nwins[i-1])
                break

        # Loop for upper bound
        for i, void in enumerate(nwins[0:-1]):
            if (nwins[i][1] < pick.fw2 and nwins[i+1][0] > pick.fw2) or \
               (nwins[i][0] <= pick.fw2 and nwins[i][1] >= pick.fw2):
                noise.append(nwins[i+1])
                break

    ampnoise = []
    for n in noise:
        startlabel = int(n[0] * 2. * np.pi/(1000. * delomeg)) + 1
        endlabel = int(n[1] * 2. * np.pi/(1000. * delomeg)) + 1
        nind.append([startlabel, endlabel])
        ampnoise.append(abs(Fxx[startlabel:endlabel+1]).max())

    ampnoise = np.array(ampnoise).max()

    startlabel = int(pick.fw1 * 2. * np.pi/(1000. * delomeg)) + 1
    endlabel = int(pick.fw2 * 2. * np.pi/(1000. * delomeg)) + 1
    ampsignal = abs(Fxx[startlabel:endlabel+1]).max()

    return ampsignal/ampnoise, noise, overlap_err


def signal2fwhm(y, x_center, max_peak='data', interpolate=True,
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
            # Check at both ends if intervals are open -> set hw to inf if so
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

    hw = abs(hw_points[1][0] - hw_points[1][-1])
    if hw == 0:
        hw = 1
        ymax = 0
    else:
        if max_peak == 'interpolated':
            ymax = env[hw_points[1]].max()
        elif max_peak == 'data':
            ymax = amp[env_shift + hw_points[1]].max()

    return ymax/hw/np.power(no_of_peaks, peak_order)


def neighbouring_minima(data, data_maxarg):
    leftfound = False
    rightfound = False
    minima = [None, None]
    i = 0
    while True:
        ri = i
        x = data[data_maxarg + i]
        xx = data[data_maxarg + ri + 1]
        if not rightfound:
            if x <= xx:
                minima[1] = int(ri)
                rightfound = True

        li = - i
        y = data[data_maxarg + li]
        yy = data[data_maxarg + li - 1]
        if not leftfound:
            if y <= yy:
                minima[0] = int(li)
                leftfound = True

        if leftfound and rightfound:
            break
        i += 1
    return minima


def takeClosest(myList, myNumber, N=1):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if N == 1:
        if pos == 0:
            return myList[0]
        if pos == len(myList):
            return myList[-1]
        before = myList[pos - 1]
        after = myList[pos]
        if after - myNumber < myNumber - before:
            return after
        else:
            return before

    if N == 2:
        try:
            before = myList[pos - 1]
            after = myList[pos]
            return before, after

        except IndexError:
            if pos == 0:
                return None, myList[0]
            if pos == len(myList):
                return None, myList[-1]
            before = myList[pos - 1]
            after = myList[pos]
            if after - myNumber < myNumber - before:
                return None, after
            else:
                return None, before


def split_digit_nondigit(x):
    dnd = []
    for f in filter(None, re.split(r'(\d+)', x)):
        dnd.append(f)
    return dnd


def mode_name(names):
    for name in names:
        n = list(name[0:5])
        # Checking for right format: Digit Letter Digit
        if n[0] == '0' and n[1].isdigit():
            n = n[1:]
            n = n.upper()

        elif n[0].isdigit() and not n[1].isdigit():
            n = n.upper()

        mode = str().join(n)
    return mode


def array2stream(ArrayData, st_original=None, network=None):
    """
    param network: Network, of with all the station information
    type network: obspy.core.inventory.network.Network
    """
    if ArrayData.ndim == 1:

        trace = obspy.core.trace.Trace(ArrayData)

        if isinstance(st_original, Stream):
            trace.stats = st_original[0].stats
        elif isinstance(st_original, Trace):
            trace.stats = st_original.stats

        return Stream(trace)

    else:
        traces = []

        for i, trace in enumerate(ArrayData):
            newtrace = obspy.core.trace.Trace(trace)
            traces.append(newtrace)

        stream = Stream(traces)

        # Just writes the network information,
        # if possible input original stream

        if isinstance(st_original, Stream):
            st_tmp = st_original.copy()
            # Checks length of ArrayData and st_original, if needed,
            # corrects trace.stats.npts value of new generated Stream-object.
            if ArrayData.shape[1] == len(st_tmp[0]):

                for i, trace in enumerate(stream):
                    trace.stats = st_tmp[i].stats

            else:

                for i, trace in enumerate(stream):
                    trace.stats = st_tmp[i].stats
                    trace.stats.npts = ArrayData.shape[1]

        elif isinstance(network, Network) and not isinstance(st_tmp, Stream):

            for trace in stream:
                trace.meta.network = network.code
                trace.meta.station = network[0].code

        return stream


def array2trace(ArrayData, st_original=None):
    if ArrayData.ndim != 1:
        try:
            stream = array2stream(ArrayData, st_original)
            return stream
        except Exception:
            msg = 'Dimension do not fit'
            raise IOError(msg)
    else:
        trace = obspy.core.trace.Trace(ArrayData)

    if isinstance(st_original, Stream):
        trace.stats = st_original[0].stats
    elif isinstance(st_original, Trace):
        trace.stats = st_original.stats

    return trace


def cat4stream(stream, client_name):

    client = Client(client_name)
    cat = obspy.core.event.Catalog()
    lat_old = None
    lon_old = None
    stime_old = None

    for i, trace in enumerate(stream):
        if hasattr(trace.stats, 'sh'):
            eventinfo = trace.stats.sh
            depth = eventinfo['DEPTH']+10
            lat = eventinfo['LAT']
            lon = eventinfo['LON']
            origin = eventinfo['ORIGIN']
            etime = origin + 300
            stime = origin - 300

        elif hasattr(trace.stats, 'ah'):
            eventinfo = trace.stats.ah.event
            depth = eventinfo.depth+10
            lat = eventinfo.latitude
            lon = eventinfo.longitude
            origin = eventinfo.origin_time
            etime = origin + 300
            stime = origin - 300

        else:
            print('trace no %i has no event information, skipped' % i)
            continue

        if i > 0 and lat == lat_old and lon == lon_old and stime == stime_old:
            continue

        cat_tmp = client.get_events(starttime=stime, endtime=etime,
                                    maxdepth=depth, latitude=lat,
                                    longitude=lon, maxradius=1)

        lat_old = lat
        lon_old = lon
        stime_old = stime

        cat.append(cat_tmp[0])

    return cat


def create_deltasignal(no_of_traces=10, len_of_traces=30000,
                       multiple=False, multipdist=2, no_of_multip=1,
                       slowness=None, zero_traces=False, no_of_zeros=0,
                       noise_level=0, non_equi=False):
    """
    function that creates a delta peak signal
    slowness = 0 corresponds to shift of 1 to each trace
    """
    if slowness:
        slowness = slowness-1
    data = np.array([noise_level * np.random.rand(len_of_traces)])

    if multiple:
        dist = multipdist
        data[0][0] = 1
        for i in range(no_of_multip):
            data[0][dist+i*dist] = 1
    else:
        data[0][0] = 1

    data_temp = data
    for i in range(no_of_traces)[1:]:
        if slowness:
            new_trace = np.roll(data_temp, slowness*i)
        else:
            new_trace = np.roll(data_temp, i)
        data = np.append(data, new_trace, axis=0)

    if zero_traces:
        first_zero = len(data)/no_of_zeros
        while first_zero <= len(data):
            data[first_zero-1] = 0
            first_zero = first_zero+len(data)/no_of_zeros

    if non_equi:
        for i in [5, 50, 120]:
            data = line_set_zero(data, i, 10)
        data, indices = extract_nonzero(data)
    else:
        indices = []

    return(data, indices)


def create_ricker(n_of_samples, n_of_traces, delta_traces=1, slope=0,
                  n_of_ricker_samples=100, width_of_ricker=2.,
                  shift_of_ricker=0):
    """
    Creates n_of_traces Traces with a Ricker wavelet
    :param n_of_samples: No of samplesw
    :type  n_of_samples: int

    :param n_of_traces: No of traces
    :type  n_of_traces: int

    :param slope: Indexshift of the traces, shift is applied by the relation
                  delta_t = delta_traces * slope
    :type  slope: int

    :param width_of_ricker: width_of_ricker parameter of Ricker-wavelet,
                            default 2
    :type  width_of_ricker: float

    :param n_of_ricker_samples: Number of samples for ricker
    :type  n_of_ricker_samples: int
    """
    if SCIPY is False:
        msg = 'Cannot find scipy library.'
        raise IOError(msg)

    if n_of_samples < n_of_ricker_samples:
        msg = 'Number of tracesamples lower than number of ricker samples'
        raise IOError(msg)

    data = np.zeros((n_of_traces, n_of_samples))

    trace = np.zeros(n_of_samples)
    ricker_tmp = sp.signal.ricker(n_of_ricker_samples, width_of_ricker)
    ricker = ricker_tmp/ricker_tmp.max()

    trace[shift_of_ricker:shift_of_ricker+n_of_ricker_samples] = ricker

    if slope != 0:
        for i in range(data.shape[0]):
            delta = np.floor(i * float(abs(slope) / float(delta_traces)))
            delta = delta.astype('int')
            data[i] = np.roll(trace, delta)[:n_of_samples]
        if slope < 0:
            data = np.flipud(data)
    elif slope == 0:
        for i, dt in enumerate(data):
            data[i] = trace

    return data


def create_sine(no_of_traces=10, len_of_traces=30000, samplingrate=30000,
                no_of_periods=1):

    deltax = 2*np.pi/len_of_traces
    signal_len = len_of_traces * no_of_periods
    data_temp = np.array([np.zeros(signal_len)])
    t = []

    # first trace
    for i in range(signal_len):
        data_temp[0][i] = np.sin(i*deltax)
        t.append((float(i) + float(i)/signal_len)*2*np.pi/signal_len)
        data = data_temp

    # other traces
    for i in range(no_of_traces)[1:]:
        data = np.append(data, data_temp, axis=0)

    return(data, t)


def cut2shortest(stream):
    """
    Cuts traces in stream to the same length. Looks for the latest beginning
    and the earliest ending of traces in stream,
    which will be the new reference times.
    """
    start = stream[0].stats.starttime
    end = stream[0].stats.endtime
    for trace in stream:
        if trace.stats.starttime > start:
            start = trace.stats.starttime
        if trace.stats.endtime < end:
            end = trace.stats.endtime

    stream.trim(start, end)
    return stream


def extract_nonzero(array):
    newarray = array[~np.all(array == 0, axis=1)]
    newindex = np.unique(array.nonzero()[0])
    return(newarray, newindex)


def fourier_transform(trace, t1, t2, n1, n2, shape, mfac=1):
    power = nextpow2(len(trace.data))
    if power < 16:
        power = 16
    maxdata = np.power(2, power)  # len(trace.data)
    indvec = np.arange(1, maxdata+1)
    delomeg = 2. * np.pi / (maxdata * trace.stats.delta)
    f = 1000. * indvec * delomeg / (2. * np.pi)
    f = f - (1000. * delomeg / (2. * np.pi))
    mid = t1+0.5*(t2-t1)
    shift_fac = np.exp(1j * f * mid * 2. * np.pi/1000.)
    tracefill = taper(n1, n2, trace, shape) * mfac

    X = np.fft.fft(tracefill, maxdata)
    FT = X * shift_fac
    return f, FT, delomeg


def taper_FT(tr, tw, fw, taper_shape='hanning'):
    wstart = fw[0]
    wend = fw[1]
    Htstart_org = tw[0] * 3600.
    Htend_org = tw[1] * 3600.
    Htstart, Htend = Htstart_org, Htend_org
    blockPrint()
    times = get_times(tr, Htstart, Htend, Htstart_org, Htend_org)
    enablePrint()
    tstart, nstart, nend, t, Htstart, Htend, Terr = times
    taper_shape = 'hanning'
    f, Fxx, delomeg = fourier_transform(tr, Htstart, Htend,
                                        nstart, nend, taper_shape)
    startlabel = int(np.round(wstart * 2. * np.pi/(1000. * delomeg)))-1
    endlabel = int(np.round(wend * 2. * np.pi/(1000. * delomeg)))-1

    return f, Fxx, startlabel, endlabel


def inv4stream(stream, client_name=None, level='channel',
               station_info='custom'):

    if client_name is not None:
        stat = []
        channels = []
        for i, trace in enumerate(stream):
            if i > 0:
                if trace.stats.station == stream[i-1].stats.station:
                    continue
            stat.append(trace.stats.station)
            channel = trace.stats.channel
            if channel.endswith('T'):
                channel = "%s%s" % (channel[:-1], 'N')
            if channel.endswith('R'):
                channel = "%s%s" % (channel[:-1], 'E')
            if station_info != 'custom':
                if channel[0] not in channels:
                    channels.append("%s*" % channel[0])
            else:
                channels.append(channel)

        stat = ','.join(list(set(stat)))
        channels = ','.join(list(set(channels)))

        client = Client(client_name)
        inv_tmp = client.get_stations(station=stat, channel=channels)
        nets = []
        for net in inv_tmp.networks:
            if net.code == 'SYS':
                continue
            nets.append(net.code)
        nets = ','.join(list(set(nets)))
        inv = client.get_stations(station=stat, network=nets,
                                  channel=channels,
                                  level=level)

    else:
        inv = Inventory(
              networks=[],
              source="nmpy")
        net = Network(
              code=".",
              stations=[])
        for tr in stream:
            if hasattr(tr.stats, 'ah'):
                station = tr.stats.ah.station
                sta = Station(
                      code=station.code,
                      latitude=station.latitude,
                      longitude=station.longitude,
                      elevation=station.elevation)
                net.stations.append(sta)
            if hasattr(tr.stats, 'sac'):
                station = tr.stats.sac
                sta = Station(
                      code=station.kstnm,
                      latitude=station.stla,
                      longitude=station.stlo,
                      elevation=station.stel)
                net.stations.append(sta)
        inv.networks.append(net)

    if len(inv) == 0:
        return None
    else:
        return inv


def inv4trace(trace, client_name=None, level='channel'):

    if client_name is not None:
        stat = str(trace.stats.station)
        channel = str(trace.stats.channel)
        if channel.endswith('T'):
            channel = "%s%s" % (channel[:-1], 'N')
        if channel.endswith('R'):
            channel = "%s%s" % (channel[:-1], 'E')
        client = Client(client_name)
        inv_tmp = client.get_stations(station=stat, channel=channel)
        nets = []
        for net in inv_tmp.networks:
            if net.code == 'SYS' or net.code == 'SY':
                continue
            nets.append(net.code)
        nets = ','.join(list(set(nets)))
        inv = client.get_stations(station=stat, network=nets,
                                  channel=channel,
                                  level=level)

    else:
        inv = Inventory(
              networks=[],
              source="nmpy")
        net = Network(
              code=".",
              stations=[])
        if hasattr(trace.stats, 'ah'):
            station = trace.stats.ah.station
            sta = Station(
                  code=station.code,
                  latitude=station.latitude,
                  longitude=station.longitude,
                  elevation=station.elevation)
            net.stations.append(sta)
        inv.networks.append(net)

    if len(inv) == 0:
        return None
    else:
        return inv


def list2stream(list):

    stream = Stream()
    for station in list:
        for trace in station:
            stream.append(trace)

    return stream


def maxrow(array):
    rowsum = 0
    for i in range(len(array)):
        if array[i].sum() > rowsum:
            rowsum = array[i].sum()
            max_row_index = i
    return(max_row_index)


def keep_longest(stream):
    """
    keeps the longest record of each channel
    """

    st_tmp = Stream()
    st_tmp.sort(['npts'])
    channels = AttribDict()

    for i, tr in enumerate(stream):

        if tr.stats.channel in channels:
            continue
        else:
            # Append the name of channel, samplingpoints and number of trace
            channels[tr.stats.channel] = [tr.stats.npts, i]
            st_tmp.append(stream[i])

    stream = st_tmp

    return stream


def nextpow2(i):
    # See Matlab documentary
    n = 1
    count = 0
    while n < abs(i):
        n = np.power(2, count)
        count += 1
    return count-1


def read_file(stream, inventory, catalog, array=False):
    """
    function to read data files, such as MSEED, station-xml and quakeml, in a
    way of obspy.read if need, pushes stream in an array for further processing
    """
    st = obspy.read(stream)
    inv = obspy.read_inventory(inventory)
    cat = obspy.readEvents(catalog)

    # pushing the trace data in an array
    if array:
        ArrayData = stream2array(st)
        return(st, inv, cat, ArrayData)
    else:
        return(st, inv, cat)


def split2stations(stream, min_len=None, merge_traces=None, keep_masked=False):
    """
    Splits a stream in a list of streams, sorted by the stations inside
    stream object. Merges traces with the same ID to one trace.

    :param stream:
    :type  stream:
    :param merge_traces: defines if traces should be merged,
                         or just the longest continious record is kept.
    :type  merge_traces: bool or none
    """
    stream.sort(['station'])

    stream_list = []
    st_tmp = Stream()

    statname = stream[0].stats.station
    for trace in stream:
        # Collect traces from same station
        if trace.stats.station == statname:
            st_tmp.append(trace)

        else:

            if merge_traces is True:
                try:
                    st_tmp.merge()
                except Exception:
                    st_tmp = keep_longest(st_tmp)
            elif merge_traces is False:
                st_tmp = keep_longest(st_tmp)

            stream_list.append(st_tmp)
            statname = trace.stats.station
            st_tmp = Stream()
            st_tmp.append(trace)

    if merge_traces is True:
        try:
            st_tmp.merge()
        except Exception:
            st_tmp = keep_longest(st_tmp)
    elif merge_traces is False:
                st_tmp = keep_longest(st_tmp)

    stream_list.append(st_tmp)

    if not keep_masked or min_len:
        for station in stream_list:
            station.sort(['channel'])
            for trace in station:
                if type(trace.data) == np.ma.core.MaskedArray:
                    stream_list.remove(station)
                    break

                elif trace.stats.npts < min_len:
                    stream_list.remove(station)
                    break

    return(stream_list)


def standard_test_signal(snes1=1, snes2=3, noise=0, nonequi=False):
    y, yindices = create_deltasignal(no_of_traces=200, len_of_traces=200,
                                     multiple=True, multipdist=5,
                                     no_of_multip=1, slowness=snes1,
                                     noise_level=noise, non_equi=nonequi)

    x, xindices = create_deltasignal(no_of_traces=200, len_of_traces=200,
                                     multiple=True, multipdist=5,
                                     no_of_multip=5, slowness=snes2,
                                     noise_level=noise, non_equi=nonequi)
    a = x + y
    y_index = np.sort(np.unique(np.append(yindices, xindices)))
    return(a, y_index)


def stats(stream):
    """
    Prints stats of the stream
    """

    for trace in stream:
        print(trace.stats)

    return


def stream2array(stream, normalize=False):
    sx = stream.copy()
    x = np.zeros((len(sx), len(sx[0].data)))
    for i, traces in enumerate(sx):
        x[i] = traces.data

    if normalize:
        if x.max() == 0:
            print('Maximum value is 0')
            return(x)

        elif math.isnan(x.max()):
            print('Maximum values are NaN, set to 0')
            n = np.isnan(x)
            x[n] = 0.

        x = x / x.max()
    return(x)


def LCM(a, b):
    """
    Calculates the least common multiple of two values
    """
    import fractions
    return abs(a * b) / fractions.gcd(a, b) if a and b else 0


def line_cut(array, shape):
    """
    Sets the array to zero, except for the 0 line and given features given
    in shape, acts as bandpass filter. "Cuts" one line out + given shape.
    For detailed information look in bowpy.filter.fk.fk_filter

    :param array: array-like
    :type  array: numpy.ndarray

    :param shape: shape and filter Information
    :type  shape: list
    """

    fil = None
    name = shape[0]
    kwarg = shape[1]
    length = shape[2]
    new_array = np.zeros(array.shape).astype('complex')
    if name in ['spike', 'Spike']:
        new_array[0] = array[0]
        return new_array

    elif name in ['boxcar', 'Boxcar'] and isinstance(length, int):
        new_array[0] = array[0]
        newrange = np.linspace(1, length, length).astype('int')
        for i in newrange:
            new_array[i] = array[i]
            new_array[new_array.shape[0]-i] = array[new_array.shape[0]-i]
        return new_array

    elif name in ['butterworth', 'Butterworth',
                  'taper', 'Taper'] and isinstance(length, int):
        fil_lh = create_filter(name, array.shape[0]/2, length, kwarg)

    elif name in ['taper', 'Taper'] and isinstance(length, int):
        fil_lh = create_filter(name, array.shape[0]/2, length, kwarg)

    fil_rh = np.flipud(fil_lh)[::-1][0:][::-1]
    fil = np.zeros(2*fil_lh.size)
    fil[:fil.size/2] = fil_lh
    fil[fil.size/2:] = fil_rh

    new_array = array.transpose() * fil
    new_array = new_array.transpose()

    return(new_array)


def line_set_zero(array, shape):
    """
    Sets line zero in array + features given in shape, acts as bandstop filter.
    For detailed information look in bowpy.filter.fk.fk_filter

    :param array: array-like
    :type  array: numpy.ndarray

    :param shape: shape and filter Information
    :type  shape: list
    """

    fil = None
    name = shape[0]
    kwarg = shape[1]
    length = shape[2]
    new_array = array

    if name in ['spike', 'Spike']:
        new_array[0] = np.zeros(array[0].size)
        return new_array

    elif name in ['boxcar', 'Boxcar'] and isinstance(length, int):
        new_array[0] = np.zeros(array[0].size)
        newrange = np.linspace(1, length, length).astype('int')
        for i in newrange:
            new_array[i] = np.zeros(array[new_array.shape[0]-i].size)
            nashape = np.zeros(array[new_array.shape[0]-i].size)
            new_array[new_array.shape[0]-i] = nashape
        return new_array

    elif name in ['butterworth', 'Butterworth',
                  'taper', 'Taper'] and isinstance(length, int):
        fil_lh = create_filter(name, array.shape[0]/2, length, kwarg)

    elif name in ['taper', 'Taper'] and isinstance(length, int):
        fil_lh = create_filter(name, array.shape[0]/2, length, kwarg)
        # fil_lh = -1. * fil_lh + 1.

    fil_rh = np.flipud(fil_lh)[::-1][1:][::-1]
    fil = np.zeros(2*fil_lh.size)
    fil[:fil.size/2] = fil_lh
    fil[fil.size/2+1:] = fil_rh
    newfil = np.ones(fil.shape)
    newfil = newfil - fil

    new_array = array.transpose() * newfil
    new_array = new_array.transpose()
    return(new_array)


def create_filter(name, length, cutoff=None, ncorner=None):

    cut = float(cutoff)/float(length)
    m = float(ncorner)

    if name in ['butterworth', 'Butterworth']:
        x = np.linspace(0, 1, length)
        y = 1. / (1. + (x/float(cut))**(2.*ncorner))

    elif name in ['taper', 'Taper']:
        shift = 0.
        fit = True
        while fit:
            cut += shift
            x = np.linspace(0, 1, length)
            y = (cut-x)*m + 0.5
            y[y > 1.] = 1.
            y[y < 0.] = 0.
            if y.max() >= 1:
                fit = False
            shift = 0.1

    else:
        msg = 'No valid name for filter found.'
        raise IOError(msg)

    return y


def _unpack_string(data):
    try:
        string = data.unpack_string().split(b'\x00',
                                            1)[0].strip().decode("utf-8")
    except Exception:
        string = data.unpack_string().split()[0].strip().decode("utf-8")
    return string


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def get_stations_from_stream(stream):
    stations = []
    for tr in stream:
        if tr.stats.station in stations:
            continue
        else:
            stations.append(tr.stats.station)

    return stations


def get_chan_from_stream(stream):
    channels = []
    for tr in stream:
        channels.append(tr.stats.channel)
    channels = list(set(channels))
    return channels


def get_net_from_stream(stream):
    nets = []
    for tr in stream:
        nets.append(tr.stats.network)
    nets = list(set(nets))
    return nets


def get_chan_from_inv(inv):
    chans = inv.get_contents()['channels']
    full_channel = list(set(chans))

    channels = []
    for c in full_channel:
        channels.append(c.split('.')[-1])

    return channels, full_channel


def stream_channels_in_inv(stream, inv):
    for tr in stream:
        try:
            station = inv.select(station=tr.stats.station)[0][0]
            station.select(channel=tr.stats.channel)[0]
        except IndexError:
            return False
    return True


def download_inv_for_stream(stream, client_name='IRIS'):

    net = get_net_from_stream(stream)
    if not net:
        net = '*'
    else:
        net = ','.join(net)

    stations = get_stations_from_stream(stream)
    stations = ','.join(stations)

    channels = get_chan_from_stream(stream)
    channels = ','.join(channels)

    client = Client(client_name)
    inventory = client.get_stations(net=net, station=stations,
                                    channel=channels, level='response')
    return inventory


def find_peaks(data, peakpick=None, closest=None):
    """
    Finds peaks in given 1D array, by search for values in data
    that are higher than the two neighbours. Except for boundary values
    :param data: 1D array-like

    :param drange: optional, range of the data distribution.

    returns:
    :param peaks: array with position on 0 axis and value on 1-axis.
    """

    pos = []
    peak = []

    # Loop through all values except first and last.
    pick = None
    for p, value in enumerate(data[1:-1]):
        if data[p+1] > data[p] and data[p+1] > data[p+2]:
            if peakpick in ['mod', 'MoD', 'Mod', 'MoP', 'Mop', 'mop']:
                if data[p] > data.mean():
                    pick = data[p]
            elif isinstance(peakpick, float) or isinstance(peakpick, int):
                if data[p] > peakpick:
                    pick = data[p]
            else:
                pick = data[p]

        if pick:  # and p != 0 and p != len(data-1):
            pos.append(p+1)
            peak.append(pick)

            pick = None

    peaks = [pos, peak]

    if closest:
        if len(peaks[1]) > 1:
            ref = closest
            diff = ref

            for i, value in enumerate(peaks[0]):
                if diff > abs(ref-value):
                    diff = abs(ref-value)
                    closest_peak = [value], [[peaks[1][i]]]

            peaks = closest_peak

    return peaks


def sort_human(l, convert_exponentials=False):
    """
    Sorts a list of strings by their actual integer value.
    e.g. a list
        ['1', '10', '2']

    will be sorted to
        ['1', '2', '10']

    param l: list to be sorted
    type  l: list of strings

    """
    def convert(text):
        return float(text) if text.isdigit() else text

    def alphanum(key):
        sk = '([-+]?[0-9]*\.?[0-9]*)'
        return [convert(c) for c in re.split(sk, key)]

    if convert_exponentials is True:
        oldnames = {}
        for l_i, x in enumerate(l):
            # Check for exponentials here
            if re.search(r'[+-]?[0-9]e[+-]', x) is not None:
                i = re.search(r'[+-]?[0-9]e[+-]', x).start()
                _istart = i
                _iend = i+3
                point_cnt = 0
                while True:
                    if not is_number(x[_istart]) and point_cnt == 1:
                        _istart = _istart + 1
                        break
                    if x[_istart] == '.':
                        point_cnt += 1
                    _istart -= 1

                while True:
                    if not is_number(x[_iend]):
                        _iend = _iend - 1
                        break
                    _iend += 1

                if float(x[_istart:_iend+1]) < 1e-10:
                    print('Sort might fail, lowest value is: 1e-10')
                new_exp = '{0:.10f}'.format(float(x[_istart:_iend+1]))
                xnew = x[:_istart] + new_exp
                xnew += x[_iend+1:]
                oldnames[xnew] = x
                l[l_i] = xnew

    # l.sort(key=alphanum)
    l = natsorted(l)

    if convert_exponentials is True:
        # replace new names with original names
        for l_i, x in enumerate(l):
            if x in oldnames:
                l[l_i] = x.replace(x, oldnames[x])

    return l


def group_consecutive_numbers(data):
    """
    example:
    data = [ 1, 4,5,6, 10, 15,16,17,18, 22, 25,26,27,28]
    cons_numbers = group_consecutive_numbers(data)

    """

    cons_numbers = []
    for k, g in groupby(enumerate(data), lambda x: x[0]-x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        cons_numbers.append((group[0], group[-1]))

    return cons_numbers


def get_times(tr, Htstart, Htend, Htstart_org, Htend_org, verbose=True):

    # Set errors to false
    STshort = False
    STlong = False
    ETlong = False
    err = False
    tstart = starttime(tr)

    # if Htstart_org is not None:
    #     Htstart, Htend = Htstart_org, Htend_org

    nstart = int(round((Htstart-tstart)/tr.stats.delta))
    nend = int(round((Htend-tstart)/tr.stats.delta)+1)
    t = np.arange(tr.stats.npts) * tr.stats.delta
    t = t + tstart

    if nstart < 0:
        Tdiff = Htend - Htstart
        Htstart = tstart
        Htend = Htstart + Tdiff
        nstart = 0
        nend = int(round((Htend-tstart)/tr.stats.delta))
        t = t + tstart
        s = Htstart / 3600.
        e = Htend / 3600.
        STshort = True

    elif nstart > len(tr.data):
        Tdiff = tr.stats.endtime - tr.stats.starttime
        Htstart = tr.stats.starttime - tr.stats.ah.event.origin_time
        Htend = Htstart + Tdiff
        nstart = 0
        nend = int(round((Htend-tstart)/tr.stats.delta))
        s = Htstart / 3600.
        e = Htend / 3600.
        STlong = True

    if nend > tr.stats.npts:
        nend = tr.stats.npts - 1
        Htend = np.floor(nend * tr.stats.delta + tstart)
        Tdiff = Htend - Htstart
        ETlong = True

    s = Htstart / 3600.
    e = Htend / 3600.

    msg = ''
    if STshort or STlong or ETlong:
        err = True
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

    return(tstart, nstart, nend, t, Htstart, Htend, err)


def starttime(trace):
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

    tstart = rsec-esec+60.*((rmn-emn) + 60.*((rhr-ehr) + 24.*(nde-ndr)))

    return tstart


def taper(nstart, nend, trace, shape='hanning'):
    npun = nend-nstart
    # ppi = 2*np.pi/(npun-1)
    # ni = npun/2
    maxdata = len(trace.data)
    if shape.lower() == 'hanning':
        # fillpart = 0.5 + 0.5 * np.cos(ppi * (np.arange(1, npun)-ni))
        fillpart = np.hanning(npun)
    elif shape.lower() == 'bartlett':
        fillpart = np.bartlett(npun)
    elif shape.lower() == 'hamming':
        fillpart = np.hamming(npun)
    elif shape.lower() == 'blackman':
        fillpart = np.blackman(npun)
    elif shape.lower() == 'kaiser':
        fillpart = np.kaiser(npun, 5)
    elif shape.lower() == 'boxcar':
        fillpart = np.ones(npun)

    fill = np.zeros(maxdata)
    fill[nstart:nend] = fillpart
    tracefill = trace.data * fill

    return tracefill


def max_sc_degrees(l):
    scdegs = np.arange(0, int(2*l)+1, 2)
    return scdegs


def max_cc_degrees(m):
    """
    m = [n1, mode1, l1, n2, mode2, l2]
    param n1, mdoe1, l1 : first mode
    param n2, mode2, l2 : second mode
    type  : list

    Usage:
    m = ['0','S','15','0','T','16']
    s = max_cc_degrees(m)

    """
    type1 = m[1]
    type2 = m[4]
    l1 = int(m[2])
    l2 = int(m[5])
    smax = abs(l1+l2)
    smin = abs(l1-l2)

    if type1.lower() == type2.lower():  # same type
        ccdegs = np.arange(smin, smax+1, 2)
    else:                              # different type
        ccdegs = np.arange(smin+1, smax, 2)

    return ccdegs


def sc_degrees(max_s):
    if max_s < 0:
        lenc = 0
    else:
        cdegs = np.arange(2, max_s+2, 2)
        lenc = 2
        for c in cdegs:
            lenc += 2*c + 1

    return lenc


def cc_degrees(m):
    """
    m = [n1, mode1, l1, n2, mode2, l2, maxse, maxsq]
    param n1, mdoe1, l1 : first mode
    param n2, mode2, l2 : second mode
    param maxse : maximum cst degree
    param maxsq : maximum dst degree
    type  : list

    Usage:
    m = ['0','S','15','0','T','16','4444','0']
    s = cc_degrees(m)

    """
    # m =
    type1 = m[1]
    type2 = m[4]
    l1 = int(m[2])
    l2 = int(m[5])
    max_ccdeg = int(m[6])
    smin = abs(l1-l2)

    if type1.lower() == type2.lower():  # same type
        ccdegs = np.arange(smin, max_ccdeg+1, 2)
    else:                              # different type
        smin = smin + 1
        ccdegs = np.arange(smin, max_ccdeg+1, 2)

    return ccdegs


def get_local_extrema(data, extrema='max'):
    """
    Function calculates local extrema of rather smooth data sets

    returns: indices of local extrema a
    """

    if extrema in ['minmax']:
        a = diff(sign(diff(data))).nonzero()[0] + 1  # local min+max
    elif extrema in ['min']:
        a = (diff(sign(diff(data))) > 0).nonzero()[0] + 1  # local min
    elif extrema in ['max']:
        a = (diff(sign(diff(data))) < 0).nonzero()[0] + 1  # local max

    return a


def mask_data(data, start, end, shape):
    L = end - start
    if shape == 'boxcar':
        data[start:end] = 0.
    elif shape == 'linear':
        x = np.arange(L, dtype=float)
        m = (data[end] - data[start]) / L
        b = data[start]
        data[start:end] = m * x + b
    elif shape == 'quad' and SCIPY is True:
        x_start = np.arange(start-2, start+1)
        y_start = data[start-2:start+1]

        x_end = np.arange(end, end+3)
        y_end = data[end:end+3]

        x = np.append(x_start, x_end)
        y = np.append(y_start, y_end)
        fxq = interp1d(x, y, kind='quadratic')

        xq = np.arange(start, end+1)
        data[start-2:end+3] = fxq(xq)
    elif shape == 'cubic' and SCIPY is True:
        x_start = np.arange(start-2, start+1)
        y_start = data[start-2:start+1]

        x_end = np.arange(end, end+3)
        y_end = data[end:end+3]

        x = np.append(x_start, x_end)
        y = np.append(y_start, y_end)
        fxc = interp1d(x, y, kind='cubic')

        xc = np.arange(start, end+1)
        data[start-2:end+3] = fxc(xc)
        # xf1 = np.arange(a.min(), a.max()+1)
    return data


def update_progress(progress, title='Percent'):
    barLength = 56  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = float(abs(progress))
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\r{0}: [{1}] {2:6.1f}% {3}"
    text = text.format(title, "#"*block + "-"*(barLength-block),
                       progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def chan_doublettes(stream):
    """
    Returns True if there is more than one channel with the same name,
            else False
    """
    channels = []

    for tr in stream:
        if tr.stats.channel in channels:
            return True
        else:
            channels.append(tr.stats.channel)
    return False


def mode_in_fw(fw1, fw2, modes):
    """
    :param fw1, fw2: Cornerpoints of frequency window in mHz
    :param modes: frospy.core.modes.Modes object
    """
    mode_list = []
    if fw1 < fw2:
        for mode in modes:
            if mode.freq <= fw2 and mode.freq >= fw1:
                mode_list.append(mode.name)
    else:
        msg = "fw2 must be bigger than fw1"
        print(msg)
    return mode_list


def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def grouper(iterable, n, fillvalue=None):
    """
    Iterate over an array (iterable) in blocks of size `n`
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def chunking_list(list, keyword):
    """
    taken from
    "http://stackoverflow.com/questions/19575702/pythonhow-to-split-
    file-into-chunks-by-the-occurrence-of-the-header-word"
    If keyword is integer, list will be chunked in sublists containing that
    many entries as defined with keyword
    e.g.
        list = ['x', 'x', 'z', 'z']
        chunking_list(list, 2) = ['x', 'x'], ['z', 'z']

    """
    chunks = []
    current_chunk = []

    if type(keyword) is int:
        for i, line in enumerate(list):
            # look for multiples of keyword
            if (i+1) % keyword == 0 and current_chunk:
                current_chunk.append(line)
                chunks.append(current_chunk[:])
                current_chunk = []
                continue
            else:
                current_chunk.append(line)
        if current_chunk:
            chunks.append(current_chunk)

    else:
        for line in list:
            if line.startswith(keyword) and current_chunk:
                chunks.append(current_chunk[:])
                current_chunk = []
            current_chunk.append(line)
        chunks.append(current_chunk)

    return chunks


def freq_index(freq, delomeg):
    return int(np.round(freq * 2. * np.pi / (1000. * delomeg))) + 1


def time_index(t, delta):
    return int(round(t/delta))


def get_nested_dict(dict, key1, key2, value):
    try:
        dict[key1][key2] = value
    except KeyError:
        dict[key1] = {key2: value}
    return dict


def merge_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


def N_splitting_coeffs(smax, smin=0):
    """
    This function calculates the number of coefficients for a given range
    of oven or odd structural degrees S (S_min to S_max).
    Default S_min is 0.

    ---------------------------------------------------------------------
    intputs:
    param smin, smax: minimum and maximum degree S
                      both even or both odd
    type smin, smax: integer

    output:
    N
    """
    # check if smin and smax are odd, for odd degree structures
    if smax < 0:
        return 0

    errmsg = "smin and smax have to be both even or odd"
    if smin % 2 == 0 and smax % 2 != 0:
        raise IOError(errmsg)
    if smin % 2 != 0 and smax % 2 == 0:
        raise IOError(errmsg)

    N = 0
    for s in np.arange(smin, smax+2, 2):
        N += (2 * s) + 1

    N = int(N)
    return N


def swap(list, i, j):
    """
    Swaps element i and j in a array-like list:

    e.g.
        x = [2, 3]
        swap(x, 0, 1): [3, 2]
    """
    list[i], list[j] = list[j], list[i]
    return list


def is_mode(string):
    x = split_digit_nondigit(string)
    if len(x) == 3:
        if (
            type(float(x[0])) == float and
            type(str(x[1])) == str and
            type(float(x[2])) == float
             ):
                if x[1].upper() == 'T':
                    return 'T'
                if x[1].upper() == 'S':
                    return 'Z'

        else:
            return False
    else:
        return False


def uniq(input):
    """
    Delete duplicate Mode in current Modes object.
    """
    output = input.copy()
    for m in input.__iter__():

        k = 0
        rem = input.__class__()
        for m2 in output.__iter__():
            if m == m2:
                k = k + 1
                rem += m2
        if k >= 2:
            for b in range(1, k):
                output.remove(rem[b])
    return output


def sort_py2(input):
    output = sorted(input, key=lambda x: (x is not None,
                                          "" if isinstance(x, Number)
                                          else type(x).__name__, x))
    return output


def uniq_modes(modes):
    output = modes.copy()
    for m in modes.__iter__():
        k = 0
        rem = modes.__class__()
        for m2 in output.__iter__():
            if m.name == m2.name:
                k = k + 1
                rem += m2
        if k >= 2:
            for b in range(1, k):
                output.remove(rem[b])
    return output


def sort_catalog(cat, value):
    cat_new = obspy.Catalog()
    if value in ('time', 'date'):
        for k in sorted(cat, key=lambda k: k.origins[0].time,
                        reverse=True):
            cat_new += k
    if value in ('longitude', 'lon'):
        for k in sorted(cat, key=lambda k: k.origins[0].logintude,
                        reverse=True):
            cat_new += k
    if value in ('latitude', 'lat'):
        for k in sorted(cat, key=lambda k: k.origins[0].latitude,
                        reverse=True):
            cat_new += k
    if value == 'depth':
        for k in sorted(cat, key=lambda k: k.origins[0].depth,
                        reverse=True):
            cat_new += k
    if value in ('magnitude', 'mag'):
        for k in sorted(cat, key=lambda k: k.magnitudes[0].mag,
                        reverse=True):
            cat_new += k
    return cat_new


def find_unique_name(filename, verbose=False):
    i = 0
    while True:
        if verbose is True:
            msg = ''
        path = os.path.dirname(filename)
        f = os.path.basename(filename)
        f, ext = os.path.splitext(f)

        if os.path.exists(filename):
            if verbose is True:
                msg += '\033[93mFile exist,'
                msg += 'not overwriting\033[0m'
            if i == 0:
                f = "%s_%s" % (f, str(i))
            i += 1
            a = "_%s" % str(i-1)
            b = "_%s" % str(i)
            f = f.replace(a, b)
            filename = os.path.join(path, "%s%s" % (f, ext))
        else:
            if verbose is True:
                print(msg)
            return filename


def fQ2cst(f, Q, mode):
    """
    f and Q in mikro Hz
    """
    f0 = mode.freq * 1e3
    Q0 = mode.Q

    Rec00 = (f - f0) * np.sqrt(4 * math.pi)
    Imc00 = np.sqrt(4. * math.pi) * (f / (2. * Q) - f0 / (2. * Q0))
    return Rec00, Imc00


def get_err(cst, zeros=False):
    if type(cst) in (float, int):
        cst = [cst]
    err = QuantityError()
    if zeros is True:
        err.uncertainty = np.zeros(len(cst))
        err.upper_uncertainty = np.zeros(len(cst))
        err.lower_uncertainty = np.zeros(len(cst))
        err.confidence_level = 0
    else:
        err.uncertainty = np.array(cst)
        err.upper_uncertainty = np.array(cst)
        err.lower_uncertainty = np.array(cst)
        err.confidence_level = 0
    return err


def fQ2cst_err(f, f_err, Q, Q_err, mode, c00, d00):
    _c00 = fQ2cst(f + f_err, Q + Q_err, mode)
    # c00 = [abs(c00[0] - cst[name]['0']), abs(c00[1] - dst[name]['0'])]
    c00_err = abs(_c00[0] - c00)
    d00_err = abs(_c00[1] - d00)
    return c00_err, d00_err


def cst2fQ(Rec00, Imc00, mode, err=None):
    """
    cst is type frospy.core.splittingfunc.splittingfunc.SplittingFunc
    """
    f0 = mode.freq * 1e3
    Q0 = mode.Q

    fc = f0 + 1. / np.sqrt(4 * math.pi) * Rec00
    if err is not None:
        Q = fc / (
            2. * (f0/(2. * (Imc00 + err) * (Q0 + 1./np.sqrt(4 * math.pi))))
                  )
    else:
        Q = fc / (
            2. * (f0 / (2. * (Q0 + 1./np.sqrt(4 * math.pi))))
                  )
    return fc, Q


def calc_Q(mode, fc, cst, err=None):
    f0 = mode.freq * 1e3
    if err is None:
        _Q = (f0 / (2*mode.Q)) + 1. / np.sqrt(4. * np.pi) * cst
    else:
        _Q = (f0 / (2*mode.Q)) + 1. / np.sqrt(4. * np.pi) * cst + err
    return 0.5 * fc / _Q


def convrate(x, L=None, i=-1):
    if L is None:
        c = abs(x[i-1] - x[i]) / abs(x[i-2] - x[i-1])
    else:
        c = abs(x[i-1] - L) / abs(x[i-2] - L)
    return c


def split2chars(word):
    return [char for char in word]
