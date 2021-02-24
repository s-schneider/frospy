from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport
from obspy.core import AttribDict
from frospy.util.base import (stream_channels_in_inv, inv4stream, _unpack_string)
from obspy import Stream, Trace, UTCDateTime

import xdrlib

import numpy as np
import os


def read_ahx(filename):
    """
    Reads an AH v1 waveform file and returns a Stream object.

    :type filename: str
    :param filename: AH v1 file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """

    def _unpack_trace(data):
        ah_stats = AttribDict({
            'version': '1.0',
            'event': AttribDict(),
            'station': AttribDict(),
            'record': AttribDict(),
            'extras': []
        })

        # station info

        ah_stats.station.code = _unpack_string(data)
        ah_stats.station.channel = _unpack_string(data)
        ah_stats.station.type = _unpack_string(data)
        ah_stats.station.latitude = data.unpack_float()
        ah_stats.station.longitude = data.unpack_float()
        ah_stats.station.elevation = data.unpack_float()
        ah_stats.station.gain = data.unpack_float()
        ah_stats.station.normalization = data.unpack_float()  # A0
        poles = []
        zeros = []
        for _i in range(0, 30):
            r = data.unpack_float()
            i = data.unpack_float()
            poles.append(complex(r, i))
            r = data.unpack_float()
            i = data.unpack_float()
            zeros.append(complex(r, i))
        # first value describes number of poles/zeros
        npoles = int(poles[0].real) + 1
        nzeros = int(zeros[0].real) + 1
        ah_stats.station.poles = poles[1:npoles]
        ah_stats.station.zeros = zeros[1:nzeros]

        # event info
        ah_stats.event.latitude = data.unpack_float()
        ah_stats.event.longitude = data.unpack_float()
        ah_stats.event.depth = data.unpack_float()
        ot_year = data.unpack_int()
        ot_mon = data.unpack_int()
        ot_day = data.unpack_int()
        ot_hour = data.unpack_int()
        ot_min = data.unpack_int()
        ot_sec = data.unpack_float()
        try:
            ot = UTCDateTime(ot_year, ot_mon, ot_day, ot_hour, ot_min, ot_sec)
        except Exception:
            ot = None
        ah_stats.event.origin_time = ot
        ah_stats.event.comment = data.unpack_string()
        # temporary fix
        ah_stats.event.comment = u''

        # record info
        ah_stats.record.type = dtype = data.unpack_int()  # data type
        ah_stats.record.ndata = ndata = data.unpack_uint()  # number of samples
        ah_stats.record.delta = data.unpack_float()  # sampling interval
        ah_stats.record.max_amplitude = data.unpack_float()
        at_year = data.unpack_int()
        at_mon = data.unpack_int()
        at_day = data.unpack_int()
        at_hour = data.unpack_int()
        at_min = data.unpack_int()
        at_sec = data.unpack_float()
        at = UTCDateTime(at_year, at_mon, at_day, at_hour, at_min, at_sec)
        ah_stats.record.start_time = at
        ah_stats.record.abscissa_min = data.unpack_float()
        ah_stats.record.comment = _unpack_string(data)
        ah_stats.record.log = _unpack_string(data)

        # extras
        ah_stats.extras = data.unpack_array(data.unpack_float)

        # unpack data using dtype from record info
        if dtype == 1:
            # float
            temp = data.unpack_farray(ndata, data.unpack_float)
        elif dtype == 6:
            # double
            temp = data.unpack_farray(ndata, data.unpack_double)
        else:
            # e.g. 3 (vector), 2 (complex), 4 (tensor)
            msg = 'Unsupported AH v1 record type %d'
            raise NotImplementedError(msg % (dtype))
        tr = Trace(np.array(temp))
        tr.stats.ah = ah_stats
        tr.stats.delta = ah_stats.record.delta
        tr.stats.starttime = ah_stats.record.start_time
        tr.stats.station = ah_stats.station.code
        tr.stats.channel = ah_stats.station.channel
        return tr

    st = Stream()
    with open(filename, "rb") as fh:
        # read with XDR library
        data = xdrlib.Unpacker(fh.read())
        # loop as long we can read records
        while True:
            try:
                tr = _unpack_trace(data)
                st.append(tr)
            except EOFError:
                break
        return st


def attach_ah_header(stream, event, inventory=None, add_response=True):
    origins = event.origins[0]
    # check if all data are in inventory:
    if inventory is None:
        inventory = inv4stream(stream, 'IRIS')
    elif not stream_channels_in_inv(stream, inventory):
        inventory = inv4stream(stream, 'IRIS')

    for tr in stream:
        station = inventory.select(station=tr.stats.station)[0][0]
        if add_response:
            channel = station.select(channel=tr.stats.channel)
        ah_stats = AttribDict({
                               'version': '1.0',
                               'event': AttribDict(),
                               'station': AttribDict(),
                               'record': AttribDict(),
                               'extras': []
                               })
        ah_stats.station.code = station.code
        ah_stats.station.channel = tr.stats.channel
        ah_stats.station.type = tr.stats.location
        ah_stats.station.latitude = station.latitude
        ah_stats.station.longitude = station.longitude
        ah_stats.station.elevation = station.elevation

        if add_response:
            resp = channel[0].response
            g = resp.instrument_sensitivity.value
            if g is None:
                ah_stats.station.gain = 0
            else:
                ah_stats.station.gain = g
            A0 = resp.response_stages[0].normalization_factor
            ah_stats.station.normalization = A0  # A0

            poles = resp.get_paz().poles
            zeros = resp.get_paz().zeros

            ah_stats.station.poles = poles
            ah_stats.station.zeros = zeros
        else:
            ah_stats.station.gain = 0
            ah_stats.station.normalization = 0
            ah_stats.station.poles = [0]
            ah_stats.station.zeros = [0]

        ah_stats.event.latitude = origins.latitude
        ah_stats.event.longitude = origins.longitude
        ah_stats.event.depth = origins.depth / 1000.
        ah_stats.event.origin_time = origins.time
        ah_stats.event.comment = ''

        ah_stats.record.type = 1  # data type
        ah_stats.record.ndata = tr.stats.npts  # number of samples
        ah_stats.record.delta = tr.stats.delta  # sampling interval
        ah_stats.record.max_amplitude = max(tr.data)
        ah_stats.record.start_time = tr.stats.starttime
        ah_stats.record.abscissa_min = 0  # Check out what this value is
        ah_stats.record.comment = ''
        ah_stats.record.log = ''
        ah_stats.extras = np.zeros(21)

        tr.stats.ah = ah_stats
    return


def write_ah(stream, filename, station=None, event=None, datatype='auto'):
    """
    Dummy function, that calls _write_ah1. Use this for the sake of a better
    function name. _write_ah1 's name comes from the obspy syntaxing.
    """
    _write_ah1(stream, filename, station, event, datatype)


def _write_ah1(stream, filename, station=None, event=None, datatype='auto'):
    """
    Writes a Stream object to an AH v1 waveform file.

    :type stream:
    :param stream: The ObsPy Stream object to write.
    :type filename: str
    :param filename: open file, or file-like object

    """
    if (
        filename.endswith('AH') or
        filename.endswith('ah') or
        filename.endswith('ahx')
       ):
        filename = os.path.splitext(filename)[0]

    def _pack_trace_with_ah_dict(tr, packer):

        # station info
        packer.pack_int(6)
        packer.pack_fstring(6, tr.stats.ah.station.code)
        packer.pack_int(6)
        try:
            packer.pack_fstring(6, tr.stats.channel)
        except Exception:
            packer.pack_fstring(6, tr.stats.ah.station.channel)

        packer.pack_int(8)
        packer.pack_fstring(8, tr.stats.ah.station.type)
        packer.pack_float(tr.stats.ah.station.latitude)
        packer.pack_float(tr.stats.ah.station.longitude)
        packer.pack_float(tr.stats.ah.station.elevation)
        packer.pack_float(tr.stats.ah.station.gain)
        packer.pack_float(tr.stats.ah.station.normalization)

        poles = tr.stats.ah.station.poles
        zeros = tr.stats.ah.station.zeros

        # Poles and Zeros
        packer.pack_float(len(poles))
        packer.pack_float(0)
        packer.pack_float(len(zeros))
        packer.pack_float(0)

        for _i in range(1, 30):
            try:
                r, i = poles[_i].real, poles[_i].imag
            except IndexError:
                r, i = 0, 0
            packer.pack_float(r)
            packer.pack_float(i)

            try:
                r, i = zeros[_i].real, zeros[_i].imag
            except IndexError:
                r, i = 0, 0
            packer.pack_float(r)
            packer.pack_float(i)

        # event info
        packer.pack_float(tr.stats.ah.event.latitude)
        packer.pack_float(tr.stats.ah.event.longitude)
        packer.pack_float(tr.stats.ah.event.depth)
        try:
            packer.pack_int(tr.stats.ah.event.origin_time.year)
            packer.pack_int(tr.stats.ah.event.origin_time.month)
            packer.pack_int(tr.stats.ah.event.origin_time.day)
            packer.pack_int(tr.stats.ah.event.origin_time.hour)
            packer.pack_int(tr.stats.ah.event.origin_time.minute)
            packer.pack_float(tr.stats.ah.event.origin_time.second)
        except Exception:
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_float(0)

        packer.pack_int(80)
        packer.pack_fstring(80, tr.stats.ah.event.comment)

        # record info
        dtype = tr.stats.ah.record.type
        packer.pack_int(dtype)
        ndata = tr.stats.npts
        packer.pack_uint(ndata)
        packer.pack_float(tr.stats.ah.record.delta)
        packer.pack_float(tr.stats.ah.record.max_amplitude)
        packer.pack_int(tr.stats.ah.record.start_time.year)
        packer.pack_int(tr.stats.ah.record.start_time.month)
        packer.pack_int(tr.stats.ah.record.start_time.day)
        packer.pack_int(tr.stats.ah.record.start_time.hour)
        packer.pack_int(tr.stats.ah.record.start_time.minute)
        packer.pack_float(tr.stats.ah.record.start_time.second)
        packer.pack_float(tr.stats.ah.record.abscissa_min)
        packer.pack_int(80)
        packer.pack_fstring(80, tr.stats.ah.record.comment)
        packer.pack_int(202)
        packer.pack_fstring(202, tr.stats.ah.record.log)

        # # extras
        packer.pack_array(tr.stats.ah.extras, packer.pack_float)

        # pack data using dtype from record info
        if dtype == 1:
            # float
            packer.pack_farray(ndata, tr.data, packer.pack_float)
        elif dtype == 6:
            # double
            packer.pack_farray(ndata, tr.data, packer.pack_double)
        else:
            # e.g. 3 (vector), 2 (complex), 4 (tensor)
            msg = 'Unsupported AH v1 record type %d'
            raise NotImplementedError(msg % (dtype))

        return packer

    def _pack_trace_wout_ah_dict(tr, packer, station, event, datatype):
        """
        Entry are packed in the same order as shown in
        _pack_trace_with_ah_dict .The missing information
        is replaced with zeros
        station info
        """
        packer.pack_int(6)
        packer.pack_fstring(6, tr.stats.station)
        packer.pack_int(6)
        packer.pack_fstring(6, tr.stats.channel)
        packer.pack_int(8)
        packer.pack_fstring(8, 'null')
        # There is no information about latitude, longitude, elevation,
        # gain and normalization in the basic stream object,  are set to 0
        if not type(station) is None and hasattr(tr.stats, 'response'):
            packer.pack_float(station.latitude)
            packer.pack_float(station.longitude)
            packer.pack_float(station.elevation)
            packer.pack_float(tr.stats.response.instrument_sensitivity.value)
            tn = tr.stats.response.response_stages[0].normalization_factor
            packer.pack_float(tn)
        else:
            packer.pack_float(0)
            packer.pack_float(0)
            packer.pack_float(0)
            packer.pack_float(0)
            packer.pack_float(0)

        # Poles and Zeros are not provided by stream object, are set to 0
        if hasattr(tr.stats, 'response'):
            poles = tr.stats.response.get_paz().poles
            zeros = tr.stats.response.get_paz().zeros
            packer.pack_float(len(poles))
            packer.pack_float(0)
            packer.pack_float(len(zeros))
            packer.pack_float(0)

            for _i in range(1, 30):
                try:
                    r, i = poles[_i].real, poles[_i].imag
                except IndexError:
                    r, i = 0, 0
                packer.pack_float(r)
                packer.pack_float(i)

                try:
                    r, i = zeros[_i].real, zeros[_i].imag
                except IndexError:
                    r, i = 0, 0
                packer.pack_float(r)
                packer.pack_float(i)
        else:
            for _i in range(0, 30):
                packer.pack_float(0)
                packer.pack_float(0)
                packer.pack_float(0)
                packer.pack_float(0)

        # event info
        if event:
            origin = event.origins[0]
            packer.pack_float(origin.latitude)
            packer.pack_float(origin.longitude)
            packer.pack_float(origin.depth)
            packer.pack_int(origin.time.year)
            packer.pack_int(origin.time.month)
            packer.pack_int(origin.time.day)
            packer.pack_int(origin.time.hour)
            packer.pack_int(origin.time.minute)
            second = origin.time.second + origin.time.microsecond / 1E6
            packer.pack_float(second)
        else:
            packer.pack_float(0)
            packer.pack_float(0)
            packer.pack_float(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_int(0)
            packer.pack_float(0)

        packer.pack_int(80)
        packer.pack_fstring(80, 'null')

        # record info
        if datatype == 'float':
            dtype = 1
        elif datatype == 'double':
            dtype = 6
        elif datatype == 'auto':
            dtype = type(tr.data[0])
            if '32' in str(dtype):
                dtype = 1
            elif '64' in str(dtype):
                dtype = 6
        packer.pack_int(dtype)
        ndata = tr.stats.npts
        packer.pack_uint(ndata)
        packer.pack_float(tr.stats.delta)
        packer.pack_float(max(tr.data))
        packer.pack_int(tr.stats.starttime.year)
        packer.pack_int(tr.stats.starttime.month)
        packer.pack_int(tr.stats.starttime.day)
        packer.pack_int(tr.stats.starttime.hour)
        packer.pack_int(tr.stats.starttime.minute)

        sec = tr.stats.starttime.second
        msec = tr.stats.starttime.microsecond
        starttime_second = float(str(sec) + '.' + str(msec))
        packer.pack_float(starttime_second)

        packer.pack_float(0)
        packer.pack_int(80)
        packer.pack_fstring(80, 'null')
        packer.pack_int(202)
        packer.pack_fstring(202, 'null')

        # # extras
        packer.pack_array(np.zeros(21).tolist(), packer.pack_float)

        # pack data using dtype from record info
        if dtype == 1:
            # float
            packer.pack_farray(ndata, tr.data, packer.pack_float)
        elif dtype == 6:
            # double
            packer.pack_farray(ndata, tr.data, packer.pack_double)
        else:
            # e.g. 3 (vector), 2 (complex), 4 (tensor)
            msg = 'Unsupported AH v1 record type %d'
            raise NotImplementedError(msg % (dtype))

        return packer

    packer = xdrlib.Packer()

    for tr in stream:
        if hasattr(tr.stats, 'ah'):
            packer = _pack_trace_with_ah_dict(tr, packer)
        else:
            packer = _pack_trace_wout_ah_dict(tr, packer, station, event,
                                              datatype)

    ofilename = filename + ".ahx"
    with open(ofilename, 'wb') as fh:
        fh.write(packer.get_buffer())
