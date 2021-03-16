from __future__ import absolute_import, print_function
from frospy.core import Spectrum
from frospy.util.array_util import (attach_network_to_traces)
from frospy.util.base import (get_stations_from_stream,
                            inv4stream, chan_doublettes, cut2shortest,
                            keep_longest)
from frospy.util.read import read_cmt_file, read_std_cat, read_st
from frospy.util.ah import attach_ah_header
import scipy
import numpy as np

import math
import sys
import os

"""
Collection of processing routines to correct data for instrument repsonse etc.
"""


def remove_response(stream, sampling_rate, pre_filt=None,
                    inventory=None, event=None, merge=False):

    # Check if response info is attached to stream
    inv_check = None
    for tr in stream:
        if not hasattr(tr.stats, 'response'):
            inv_check = inventory

        if tr.stats._format == 'AH':
            tr.stats.location = tr.stats.ah.station.type

    inventory = inv_check

    stream.detrend('constant')
    stream.detrend('linear')
    if not pre_filt:
        pre_filt = [0.00015, 0.0002, .5, .51]

    for tr in stream:
        try:
            if inventory:
                attach_network_to_traces(stream, inventory)
                tr.remove_response(inventory=inventory, pre_filt=pre_filt,
                                   output="ACC", water_level=60)
                if tr.stats.sampling_rate != sampling_rate:
                    tr.resample(sampling_rate)
            else:
                tr.remove_response(pre_filt=pre_filt,
                                   output="ACC", water_level=60)
        except ValueError:
            print('\n\ncheck_channel: Illegal RESP format\n')
            print('removed station: \n')
            stw = stream.select(station=tr.stats.station)
            print(stw)
            print('\n')

            for t in stw:
                stream.remove(t)

            msg = "Error message: \n%s\n%s\n" % (tr, sys.exc_info()[1])
            with open('error.log', 'a') as logfile:
                logfile.write(msg)
        except Exception:
            print("\n\nERROR in RESP")
            stw = stream.select(station=tr.stats.station)
            print('\n')

            for t in stw:
                stream.remove(t)

            msg = "Error message: \n%s\n%s\n" % (tr, sys.exc_info()[1])
            with open('error.log', 'a') as logfile:
                logfile.write(msg)

    stream.detrend('constant')
    stream.detrend('linear')
    if merge:
        stream.merge()

        for tr in stream:
            if type(tr.data) == np.ma.core.MaskedArray:
                stream.remove(tr)

    return stream


def rotate(stream_input, event=None, inv=None, kind='NE2RT',
           verbose=True):
    """
    Rotates traces in stream.
    kind: 'NE2RT', "ABC2T"
    """
    stream = stream_input.copy()
    stations = get_stations_from_stream(stream)

    for station in stations:
        st_tmp = stream.select(station=station)
        if kind == 'NE2RT' and len(st_tmp) < 2:
            if verbose:
                msg = "Too few components, removing:"
                print(msg)
                print(st_tmp)
            stream.remove(st_tmp[0])
            continue

        elif hasattr(st_tmp[0].stats, 'ah'):
            if event:
                src_lat = event.origins[0].latitude
                src_lon = event.origins[0].longitude
            else:
                src_lat = st_tmp[0].stats.ah.event.latitude
                src_lon = st_tmp[0].stats.ah.event.longitude
            rcv_lat = st_tmp[0].stats.ah.station.latitude
            rcv_lon = st_tmp[0].stats.ah.station.longitude

        else:
            try:
                net = st_tmp[0].stats.network
                name = st_tmp[0].stats.station
                chan = st_tmp[0].stats.channel
                s = "%s.%s..%s" % (net, name, chan)
                t = st_tmp[0].stats.starttime
                c = inv.get_coordinates(s, t)

                rcv_lat = c['latitude']
                rcv_lon = c['longitude']
                src_lat = event.origins[0].latitude
                src_lon = event.origins[0].longitude
            except Exception:
                # msg = 'No coordinates for %s' % st_tmp[0]
                continue

        # Geocentric correction is performed in get_baz_synseis
        # if correction == 'delaz':
        #     delta, azep, azst = geo.delaz(src_lat, src_lon, rcv_lat, rcv_lon)
        #     baz = azep

        baz = get_baz_synseis(src_lat, src_lon, rcv_lat, rcv_lon)

        if kind.upper() in ['NE2RT', 'RT2NE']:
            # r, t = rotate_ne_rt(n, e, baz)
            if kind.upper() == 'RT2NE':
                tr1 = st_tmp.select(channel='*R*')[0]
                tr2 = st_tmp.select(channel='*T*')[0]
                r = tr1.data
                t = tr2.data
                n, e = rotate_ne2rt(r, t, -baz)
                tr1.data = n
                tr2.data = e
                tr1.stats.channel = tr1.stats.channel.replace('R', 'N')
                tr2.stats.channel = tr2.stats.channel.replace('T', 'E')
            else:
                tr1 = st_tmp.select(channel='*N*')[0]
                tr2 = st_tmp.select(channel='*E*')[0]
                n = tr1.data
                e = tr2.data
                r, t = rotate_ne2rt(n, e, baz)
                tr1.data = r
                tr2.data = t
                tr1.stats.channel = tr1.stats.channel.replace('N', 'R')
                tr2.stats.channel = tr2.stats.channel.replace('E', 'T')

            if hasattr(st_tmp[0].stats, 'ah'):
                for tr in st_tmp:
                    channel = tr.stats.channel
                    tr.stats.ah.station.channel = channel
        elif kind.upper() == 'ABC2T':
            if inv is None:
                raise IOError('Inventory must be given')
            try:
                A = get_ABC_trace(st_tmp, '*A')[0].data
                alphaA = get_ABC_azi(inv, '*A')

                B = get_ABC_trace(st_tmp, '*B')[0].data
                alphaB = get_ABC_azi(inv, '*B')

                C = get_ABC_trace(st_tmp, '*C')[0].data
                alphaC = get_ABC_azi(inv, '*C')
            except IndexError:
                A = get_ABC_trace(st_tmp, '*1')[0].data
                alphaA = get_ABC_azi(inv, '*1')

                B = get_ABC_trace(st_tmp, '*2')[0].data
                alphaB = get_ABC_azi(inv, '*2')

                C = get_ABC_trace(st_tmp, '*3')[0].data
                alphaC = get_ABC_azi(inv, '*3')

            b = np.array([A,B,C])
            Amat = np.array([
            [np.cos(alphaA)**2, np.sin(alphaA)**2, -np.sin(2*alphaA)],
            [np.cos(alphaB)**2, np.sin(alphaB)**2, -np.sin(2*alphaB)],
            [np.cos(alphaC)**2, np.sin(alphaC)**2, -np.sin(2*alphaC)]
            ])
            xx, yy, xy = scipy.linalg.solve(Amat, b)


            Tmat = np.array([
            [np.cos(alphaA + (baz-90.-alphaA))**2, np.sin(alphaA + (baz-90.-alphaA))**2, -np.sin(2*alphaA + (baz-90.-alphaA))],
            [np.cos(alphaB + (baz-90.-alphaB))**2, np.sin(alphaB + (baz-90.-alphaB))**2, -np.sin(2*alphaB + (baz-90.-alphaB))],
            [np.cos(alphaC + (baz-90.-alphaC))**2, np.sin(alphaC + (baz-90.-alphaC))**2, -np.sin(2*alphaC + (baz-90.-alphaC))]
            ])

            Tsignal = np.zeros(len(xx))

            for row in Tmat:
                Tsignal += row[0] * xx + row[1] * yy + row[2] * xy

            Tsignal = Tsignal / 3.

            tr_new = st_tmp[0].copy()
            tr_new.data = Tsignal
            tr_new.stats.channel = 'T'
            stream += tr_new
        # print('Station: %s' % station)
        # print('src coordinates %f %f' % (src_lat, src_lon))
        # print('rcv coordinates %f %f' % (rcv_lat, rcv_lon))
        # print('BAZ %f' % baz)

    return stream


def get_ABC_trace(stream, channel):
    trace = stream.select(component=channel)
    if len(trace) > 1:
        print(trace)
        raise IOError('More than 1 A channel found')

    return trace


def get_ABC_azi(inv, channel):
    station = inv.networks[0].select(channel=channel)
    for _s in station:
        for _i, chan in enumerate(_s.channels):
            if _i == 0:
                azi = chan.azimuth
            else:
                if azi != chan.azimuth:
                    print(_s.channels)
                    raise IOError('Channels differ in azimuth!')
    return azi


def rotate_cmt(stream_input, cmt_file, correction=None, verbose=True):
    """
    Rotates traces in stream.
    """
    stream = stream_input.copy()
    stations = get_stations_from_stream(stream)
    cmt = read_cmt_file(cmt_file, full=True)

    for i, station in enumerate(stations):
        st_tmp = stream.select(station=station)
        if len(st_tmp) < 2:
            if verbose:
                msg = "Too few components, removing:"
                print(msg)
                print(st_tmp)
            stream.remove(st_tmp[0])
            continue

        src_lat = cmt[0][0]
        src_lon = cmt[0][1]
        rcv_lat = st_tmp[0].stats.ah.station.latitude
        rcv_lon = st_tmp[0].stats.ah.station.longitude

        baz = get_baz_synseis(src_lat, src_lon, rcv_lat, rcv_lon)

        tr1 = st_tmp.select(channel='*N*')[0]
        tr2 = st_tmp.select(channel='*E*')[0]
        n = tr1.data
        e = tr2.data

        r, t = rotate_ne2rt(n, e, baz)

        tr1.data = r
        tr2.data = t

        tr1.stats.channel = tr1.stats.channel.replace('N', 'R')
        tr2.stats.channel = tr2.stats.channel.replace('E', 'T')

        if hasattr(st_tmp[0].stats, 'ah'):
            for tr in st_tmp:
                channel = tr.stats.channel
                tr.stats.ah.station.channel = channel

        # print('Station: %s' % station)
        # print('src coordinates %f %f' % (src_lat, src_lon))
        # print('rcv coordinates %f %f' % (rcv_lat, rcv_lon))
        # print('BAZ %f' % baz)

    return stream


def check_delta_start_origin(stream, event, delta_h):
    """
    Adjusts the lengths of each trace to the shortest trace of that station
    If event is given, data that startes more than 1 hour after origin time
    will be deleted
    """

    stations = get_stations_from_stream(stream)

    if event:
        ot = event.origins[0].time
    else:
        ot = None

    for station in stations:
        st_tmp = stream.select(station=station)
        if ot is not None:
            for tr in st_tmp:
                tdiff = abs(ot - tr.stats.starttime)
                if tdiff > delta_h * 3600.:
                    stream.remove(tr)

    return stream


def check_stream_stat_doublettes(st):
    stations = get_stations_from_stream(st)
    for s in stations:
        stt = st.select(station=s)
        if chan_doublettes(stt):
            msg = "\n\nFound %i component(s):\n" % len(stt)
            for tr in stt:
                msg += "%s\n" % tr.__original_str__()
            with open('error.log', 'a') as logfile:
                logfile.write(msg)
            st = keep_longest(stt)
    return st


def get_baz_synseis(src_lat, src_lon, rcv_lat, rcv_lon, correction=False):
    """
    Calculating the baz according to lines 704 - 738 in
    subroutine compute_synt in
    Codes/synseis_hybrid/recv_terms.f90

    """
    lat1 = geocentric(src_lat)
    lon1 = src_lon
    lat2 = geocentric(rcv_lat)
    lon2 = rcv_lon
    lat1 = 90. - lat1  # Colatitude
    lat2 = 90. - lat2  # Colatitude

    # Calculate bazimuth using eulers formula
    lat1 = np.float64(lat1/180.*np.pi)
    lon1 = np.float64(lon1/180.*np.pi)
    lat2 = np.float64(lat2/180.*np.pi)
    lon2 = np.float64(lon2/180.*np.pi)

    print("NOTE: WE HAVE TO CHECK EULER, IT IS UPDATED FOR PYTHON 3")
    ba_rad = euler(lat1, lon1, lat2, lon2)

    beta = -ba_rad[0] * 180. / np.pi
    baz = np.mod(360. - beta, 360.)
    return baz


def rotate_ne2rt(north, east, baz):
    """
    Calculating the rotation of data according to lines 250 - 256 in
    subroutine reprema in
    Codes/synseis_hybrid/recv_terms.f90

    """
    # north = -north
    R = -north*np.cos((baz)*2.*np.pi/360.) - east * np.sin((baz)*2.*np.pi/360.)
    T = north*np.sin((baz)*2.*np.pi/360.) - east * np.cos((baz)*2.*np.pi/360.)

    return R, T


def geocentric(lat):
    # Correction according to line 2545 in compu-deriv_inv.f
    geoco=0.993277
    rad=2.0*np.pi/360.0
    lat_corr = np.arctan(geoco*np.tan(lat*rad))/rad

    return lat_corr


def euler(t1, p1, t2, p2):

    prod1 = np.zeros((3, 3))
    prod2 = np.zeros((3, 3))
    rot = np.zeros((3, 3))

    rt1 = np.array([[math.cos(t1), 0, math.sin(t1)],
                    [0, 1, 0],
                    [-math.sin(t1), 0, math.cos(t1)]
                    ])

    rt2 = np.array([[math.cos(t2), 0, -math.sin(t2)],
                    [0, 1, 0],
                    [math.sin(t2), 0, math.cos(t2)],
                    ])

    rp1 = np.array([[math.cos(p1), -math.sin(p1), 0],
                    [math.sin(p1), math.cos(p1), 0],
                    [0, 0, 1]
                    ])

    rp2 = np.array([[math.cos(p2), math.sin(p2), 0],
                    [-math.sin(p2), math.cos(p2), 0],
                    [0, 0, 1]
                    ])
    # rt1 = rt1.transpose()
    # rt2 = rt2.transpose()
    # rp1 = rp1.transpose()
    # rp2 = rp2.transpose()

    for i in range(3):
        for j in range(3):
            for k in range(3):
                prod1[i, j] = prod1[i, j] + rp2[i, k] * rt2[k, j]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                prod2[i, j] = prod2[i, j] + rp1[i, k] * prod1[k, j]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                rot[i, j] = rot[i, j]+rt1[i, k]*prod2[k, j]

    cb = rot[2, 2]
    if (cb >= 1.):
        b = 0.
    elif (cb <= -1.0):
        b = np.pi
    else:
        b = np.arccos(cb)

    sb = np.sin(b)
    if (abs(sb) <= 1e-15):
        a = p2-p1
        g = 0.
    else:
        ca = rot[2, 0] / sb
        sa = rot[2, 1] / sb
        if (abs(ca-1.0) <= 1e-8):
            a = 0.0
        else:
            if (abs(ca+1.0) <= 1e-8):
                a = np.pi
            else:
                a = np.arccos(ca)

        if (sa < (0.0)):
            a = -1.0 * a

        cg = -rot[0, 2] / sb
        sg = rot[1, 2] / sb

        if (abs(cg-1.0) < 1e-8):
            g = 0.0
        else:
            if (abs(cg+1.0) < 1e-8):
                g = np.pi
            else:
                g = np.arccos(cg)

        if (sg < (0.0)):
            g = -1.0 * g

    return a, b, g


def correct_stream(st, event, inv='IRIS', cmt_id=None, ofile=None,
                   channel='VH',
                   keep_longest_traces=False,
                   cut_to_same_length=False,
                   rotate_traces=False,
                   rm_response=False,
                   rm_tidal=False,
                   rm_tidal_path='/net/home/deuss/bin/remtidah'):

    """
    Documentaion follows
    """

    pmsg = 'Processing Data'
    print('%s' % pmsg, end='\r')
    sys.stdout.flush()

    if inv == 'IRIS':
            respmsg = pmsg + ': Fetching response information from IRIS ...'
            print('%s' % respmsg, end='\r')
            sys.stdout.flush()
            inv = inv4stream(st, 'IRIS')
            respmsg = pmsg + ':                                            '
            print('%s' % respmsg, end='\r')
            sys.stdout.flush()
    if type(event) == str:
        cmt_id = event
        event = read_std_cat(cmt_id)[0]

    if channel == 'VH':
        sampling_rate = 0.1

    if cmt_id is None:
        cmt_id = 'corrected_stream'

    if keep_longest_traces:
        st = keep_longest(st)
    if cut_to_same_length:
        st = cut2shortest(st)

    if rotate_traces:
        rotmsg = pmsg + ': Rotating Data ...'
        print('%s' % rotmsg, end='\r')
        sys.stdout.flush()
        try:
            cmt_path = '//nfs/stig/simons/alldata/cmts'
            cmt_file = "%s/%s.cmt" % (cmt_path, cmt_id)
            event = read_std_cat(cmt_file)[0]
            attach_ah_header(st, event, inv, add_response=False)
            st = rotate_cmt(st, cmt_file=cmt_file,
                            correction='geocentric')
            rotmsg1 = rotmsg + ' Using local cmt solution'
            print('%s' % rotmsg1, end='\r')
            sys.stdout.flush()
        except Exception:
            event = event
            st = rotate(st, inv=inv, event=event,
                        correction='geocentric')
            rotmsg2 = rotmsg + ' Using globalcmt solution'
            print('%s' % rotmsg2, end='\r')
            sys.stdout.flush()

    if rm_response:
        rmsg = pmsg + ': Removing response ...'
        print('%s' % rmsg, end='\r')
        sys.stdout.flush()
        st = remove_response(st, sampling_rate=sampling_rate, event=event,
                             pre_filt=[0.00015, 0.0002, .5, .51])

        st.sort(['channel'])
        st.sort(['station'])

    if rm_tidal:
        rmtidemsg = pmsg + ': Removing tidal signal...'
        print('%s' % rmtidemsg, end='\r')
        sys.stdout.flush()

        ofile_tmp1 = '%s.wtide.ahx' % cmt_id
        ofile_tmp2 = '%s.wotide.ahx' % cmt_id
        st.write(ofile_tmp1, format='AH')
        # Using the fortran code defined in rm_tidal_path
        os.system('%s %s %s' % (rm_tidal_path, ofile_tmp1, ofile_tmp2))
        st = read_st(ofile_tmp2, format='AH')
        os.system('rm %s %s' % (ofile_tmp1, ofile_tmp2))

    return st


def select_stations(infile, tw=[5, 60], min_snr=1.2,
                    verbose=False, summary_file=False):
    """
    Loops through all stations in infile and calulcates a mean snr for the
    frequencies given in fw in 1mHz intervals.

    """
    def get_snr_between_nwins(Fxx, swin, nw1, nw2):
        signal = Fxx[spec.flabel(swin[0]):spec.flabel(swin[1])+1]
        signal = signal.max()
        noise = Fxx[spec.flabel(nw1[0]):spec.flabel(nw1[1])+1]
        noise = noise.max()
        noise2 = Fxx[spec.flabel(nw2[0]):spec.flabel(nw2[1])+1]
        noise2 = noise2.max()

        if noise2 > noise:
            noise = noise2
        return signal/noise

    from nmpy import data as nmpydata
    path = nmpydata.__path__[0] + "/AD/noisewindows_new.dat"
    nwins = np.genfromtxt(path)
    st = read_st(infile, 'ah')
    st.sort()

    if summary_file is not False:
        cmt = infile.split('/')[-1].split('.ahx')[0]
        msg = '# Signal to Noise Ratios for each frequency window\n\n'
        msg += '%s, 0-1mHz, 1-2mHz, 2-3mHz, 3-4mHz, 4-5mHz, ' % cmt
        msg += '5-6mHz, 6-7mHz, 7-8mHz, 8-9mHz, 8-10mHz\n'
    for tr in st:
        spec = Spectrum(tr, tw[0], tw[1])
        Fxx = abs(list(spec.data.fft.values())[0])

        npairs = zip(nwins, nwins[1:])
        snr = []
        keep = False

        snr_dict = {}

        for nw1, nw2 in npairs:
            swin = [nw1[1], nw2[0]]
            snr = get_snr_between_nwins(Fxx, swin, nw1, nw2)
            fwin = int(np.ceil(swin[0]))

            if fwin in snr_dict:
                if snr > snr_dict[fwin]:
                    snr_dict[fwin] = snr
            else:
                snr_dict[fwin] = snr

            if snr >= min_snr:
                keep = True

        if summary_file is not False:
            if keep is True:
                msg += '%s, ' % tr.stats.station
                for i in range(1, 11):
                    if i in snr_dict:
                        msg += '%2.2f, ' % snr_dict[i]
                    else:
                        msg += '-, '
                msg += '\n'

        if keep is False:
            if verbose is True:
                print('%s removed' % tr.stats.station)
            st.remove(tr)

    if summary_file is not False:
        if verbose:
            print(msg)
        with open(summary_file, 'w') as fh:
            fh.write(msg)

    return st
