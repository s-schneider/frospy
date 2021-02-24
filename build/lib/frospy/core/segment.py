# -*- coding: utf-8 -*-
"""
Module for handling nmPy segment objects.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQAcd
import glob
import fnmatch
import os
import sys
import re
import json
from copy import deepcopy
from frospy import data as frospydata

from frospy.util.write import write_pickle
from frospy.util.read import read_pickle

from frospy.core.pick import Pick
from frospy.core.modes import read as read_modes
from frospy.core.database.query import segments_query
from frospy.util.base import swap, find_unique_name

import sqlite3


def convert_to_dat(filename):
    segment = read(filename)
    cmt_id = segment[0].stats.event
    segment.write("%s.dat", cmt_id)
    return


def read(ifile, format=None, **kwargs):
    """
    ifile: filenam, AD_picks or db
    """
    if 'modes' in kwargs:
        if type(kwargs['modes']) != list:
            modes = [kwargs['modes']]
        else:
            modes = kwargs['modes']
    else:
        modes = 'all'

    accepted_cmt = ['010196C',
                    '010498A',
                    '011301C',
                    '011307A',
                    '011593C',
                    '012601A',
                    '012706C',
                    '021796B',
                    '022710A',
                    '022710Q',
                    '030385A',
                    '030390B',
                    '030477A',
                    '030684B',
                    '030688B',
                    '030994E',
                    '031111B',
                    '031111Q',
                    '031883A',
                    '032598B',
                    '032800C',
                    '032805D',
                    '032805Q',
                    '040107E',
                    '041890B',
                    '042006K',
                    '042291A',
                    '050306F',
                    '050515A',
                    '050786B',
                    '051208A',
                    '051208B',
                    '052389A',
                    '052413A',
                    '052581A',
                    '052683A',
                    '053015A',
                    '060294C',
                    '060400D',
                    '060994A',
                    '060994Q',
                    '061096B',
                    '061278A',
                    '061796A',
                    '061796Q',
                    '061800A',
                    '062277A',
                    '062282A',
                    '062301E',
                    '062301Q',
                    '070508A',
                    '070508Q',
                    '070701F',
                    '071293B',
                    '071503F',
                    '071690A',
                    '071706B',
                    '071780A',
                    '072916A',
                    '073095A',
                    '080403C',
                    '080807G',
                    '081412A',
                    '081507F',
                    '081676B',
                    '081799A',
                    '081902C',
                    '081918A',
                    '081977B',
                    '082712A',
                    '083112A',
                    '090292A',
                    '090512A',
                    '090802H',
                    '090817A',
                    '090905C',
                    '091615A',
                    '092099D',
                    '092503C',
                    '092605B',
                    '092807I',
                    '100483C',
                    '100494B',
                    '100494Q',
                    '100805A',
                    '100995C',
                    '101192C',
                    '101497A',
                    '102086B',
                    '102615A',
                    '110302J',
                    '110897A',
                    '111296D',
                    '111401B',
                    '111407B',
                    '111506F',
                    '111703B',
                    '111713A',
                    '112084A',
                    '112483A',
                    '112978A',
                    '112998B',
                    '113076A',
                    '113087D',
                    '120395E',
                    '120597C',
                    '120678A',
                    '120907C',
                    '121279A',
                    '121292B',
                    '122291B',
                    '122304A',
                    '122516A',
                    '122604A',
                    '122604Q',
                    '122894C',
                    '123090D']

    if ifile == 'AD_picks':
        f = frospydata.__path__[0] + "/AD/AD_picks.json"
        with open(f) as fh:
            data = json.load(fh)

        all_modes = read_modes()
        segments = Segment()

        for mode, values in data['segments'].items():
            if modes != 'all':
                if mode not in modes:
                    continue

            m = all_modes.select(name=mode)
            for cmt, picks in values.items():
                for chan, pdata in picks.items():
                    for stat,  p in pdata.items():
                        fw1 = p['fw1']
                        fw2 = p['fw2']
                        tw1 = p['tw1']
                        tw2 = p['tw2']
                        weight = p['weight']
                        pick = Pick({'station': stat,
                                     'event': cmt,
                                     'channel': chan,
                                     'modes': m,
                                     'author': 'Arwen Deuss',
                                     'data_origin': 'AD_picks.json'},
                                    fw1, fw2, tw1, tw2, weight)
                        # if mode == '0S6':
                        #     print(pick.stats.modes)
                        segments += pick

        return segments

    elif ifile in ('Z', 'T', 'R'):
        f = frospydata.__path__[0] + "/SAS/VH{c}_picks.pickle".format(c=ifile)
        return read_pickle(f)

    elif ifile in ('db', 'database'):
        if 'path' in kwargs:
            f = kwargs['path']
            kwargs.pop('path')
        else:
            f = frospydata.__path__[0] + "/SAS/picks.sqlite3"

        channel = kwargs['channel']
        modes = kwargs['modes']
        min_snr = None
        author = None
        event = None
        station = None

        if 'min_snr' in kwargs:
            min_snr = kwargs['min_snr']
        if 'author' in kwargs:
            author = kwargs['author']
        if 'event' in kwargs:
            event = kwargs['event']
        if 'station' in kwargs:
            station = kwargs['station']

        dbq = segments_query(f, channel=channel, modes=modes, author=author,
                             min_snr=min_snr, event=event, station=station)
        segments = Segment()

        event_stat = {}

        for i, p in enumerate(dbq):
            name = p[0]
            fw1 = p[1]
            fw2 = p[2]
            tw1 = p[3]
            tw2 = p[4]
            weight = p[5]
            event = p[6]
            modes = p[7]
            snr = p[8]
            author = p[9]
            data_origin = p[10]

            if event not in accepted_cmt:
                if 'events' in kwargs:
                    if event not in kwargs['events']:
                        continue
                else:
                    continue

            if event not in event_stat:
                event_stat[event] = []
            if name in event_stat[event]:
                continue

            pick = Pick({
                        'station': name, 'event': event, 'author': author,
                        'data_origin': data_origin, 'modes': [modes],
                        'snr': snr, 'channel': channel},
                        fw1, fw2, tw1, tw2, weight)

            event_stat[event].append(name)

            segments += pick
        segments.del_duplicates()
        segments.sort()

        return segments

    if format is None:
        # try to guess format from file extension
        _, format = os.path.splitext(ifile)
        format = format[1:]

    if format == 'pickle':
        files = glob.glob(ifile)
        segments = Segment()
        for f in files:
            x = read_pickle(f)
            segments += update_segment_pickle(x)
        return segments

    else:
        files = glob.glob(ifile)
        segments = Segment()
        for f in files:
            with open(f, 'r') as fh:
                segment_file = fh.readlines()

            # Try to guess event id from segment file name
            segfile = f.split('/')[-1]
            try:
                match = re.search(r"[0-9]{6}[A-Z]{1}", segfile)
                evname = match.group(0)
            except AttributeError:
                print('no event name assumed')

            # Create empty list with information of each station
            station = []

            for _i, line in enumerate(segment_file):
                station.append(line.split())
                # Read in blocks of 4 lines
                count = _i + 1
                if _i > 0 and count % 4 == 0:
                    name = station[0][0]
                    fw1 = float(station[1][0])
                    fw2 = float(station[1][1])
                    tw1 = float(station[2][0])
                    tw2 = float(station[2][1])
                    weight = float(station[3][0])
                    pick = Pick({'station': name, 'event': evname}, fw1, fw2,
                                tw1, tw2, weight)
                    segments += pick
                    station = []

        return segments


def update_segment_pickle(segment):
    segment_new = Segment()

    for key, value in vars(segment).items():
        setattr(segment_new, key, value)

    return segment_new


class Segment(object):
    """
    Container for frospy.core.pick objects
    """
    def __init__(self, picks=None):
        self.picks = []
        if isinstance(picks, Pick):
            picks = [picks]
        if picks is not None:
            self.picks.extend(picks)

    def __str__(self, extended=False):
        out = str(len(self.picks)) + ' Picks(s) in Segment:\n'
        if len(self.picks) <= 20 or extended is True:
            out = out + "\n".join([_i.__str__() for _i in self])
        else:
            out = out + "\n" + self.picks[0].__str__() + "\n" + \
                '...\n(%i other picks)\n...\n' % (len(self.picks) - 2) + \
                self.picks[-1].__str__() + '\n\n[Use "print(' + \
                'Segment.__str__(extended=True))" to print all Traces]'
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def __nonzero__(self):
        """
        A Segment is considered zero if has no Picks.
        """
        return bool(len(self.picks))

    def __len__(self):
        """
        Return the number of Picks in the Segment object.

        .. rubric:: Example

        >>> segment = Segment([Pick(), Pick(), Pick()])
        >>> len(segment)
        3
        """
        return len(self.picks)

    def __iter__(self):
        """
        Return a robust iterator for segment.picks.

        Doing this it is safe to remove picks from Segments inside of
        for-loops using segment's :meth:`~frospy.core.segment.Segment.remove`
        method. Actually this creates a new iterator every time a pick is
        removed inside the for-loop.

        """
        return list(self.picks).__iter__()

    def __add__(self, other):
        """
        Add two Segment or a Segment with a single Pick.

        :type other: :class:`~frospy.core.segment.Segment` or
            :class:`~frospy.core.pick.Pick`
        :param other: Segment or Pick object to add.
        :rtype: :class:`~frospy.core.segments.Segment`
        :returns: New Segment object containing references to the Picks of the
            original Segment

        .. rubric:: Examples

        1. Adding two Segment

            >>> st1 = Segment([Pick(), Pick(), Pick()])
            >>> len(st1)
            3
            >>> st2 = Segment([Pick(), Pick()])
            >>> len(st2)
            2
            >>> Segment = st1 + st2
            >>> len(Segment)
            5

        2. Adding Segment and Pick

            >>> Segment2 = st1 + Pick()
            >>> len(Segment2)
            4
        """
        if isinstance(other, Pick):
            other = Segment([other])
        if not isinstance(other, Segment):
            raise TypeError
        picks = self.picks + other.picks
        return self.__class__(picks=picks)

    def __iadd__(self, other):
        """
        Add two segments with self += other.

        It will extend the current Stream object with the traces of the given
        Stream. Traces will not be copied but references to the original traces
        will be appended.

        :type other: :class:`~frospy.core.segment.Segment` or
            :class:`~obspy.core.pick.Pick`
        :param other: Segment or Pick object to add.

        .. rubric:: Example

        >>> segment = Segment([Pick(), Pick(), Pick()])
        >>> len(segment)
        3

        >>> segment += Segment([Pick(), Pick()])
        >>> len(segment)
        5

        >>> segment += Pick()
        >>> len(segment)
        6
        """
        if isinstance(other, Pick):
            other = Segment([other])
        if not isinstance(other, Segment):
            raise TypeError
        self.extend(other.picks)
        return self

    def __getitem__(self, index):
        """
        __getitem__ method of frospy.Segment objects.

        :return: Pick objects
        """
        return self.picks.__getitem__(index)

    def append(self, pick):
        """
        Append a single Pick object to the current Segment object.
        """
        if isinstance(pick, Pick):
            self.picks.append(pick)
        else:
            msg = 'Append only supports a single pick object as an argument.'
            raise TypeError(msg)
        return self

    def copy(self):
        return deepcopy(self)

    def del_duplicates(self, exact=True):
        """
        Delete duplicate picks in current Segment object.
        Duplicate picks have the same time window, frequency window,
        weight and station.
        """
        cop = self.copy()
        for pick1 in self.__iter__():
            k = 0
            rem = Segment()
            for pick2 in cop.__iter__():
                if exact is True:
                    if pick1 == pick2:
                        k = k + 1
                        rem.append(pick2)
                else:
                    if (
                        pick1.stats.event == pick2.stats.event and
                        pick1.stats.station == pick2.stats.station
                       ):
                        k = k + 1
                        rem.append(pick2)
            if k >= 2:
                for b in range(1, k):
                    cop.remove(rem[b])

        self = cop
        return

    def extend(self, pick_list):
        """
        Extend the current Segment object with a list of Pick objects.
        """
        if isinstance(pick_list, list):
            for _i in pick_list:
                # Make sure each item in the list is a trace.
                if not isinstance(_i, Pick):
                    msg = 'Extend only accepts a list of Pick objects.'
                    raise TypeError(msg)
            self.picks.extend(pick_list)
        elif isinstance(pick_list, Segment):
            self.picks.extend(pick_list.picks)
        else:
            msg = 'Extend only supports a list of Pick objects as argument.'
            raise TypeError(msg)
        return self

    def remove(self, pick):
        """
        Remove the first occurrence of the specified Pick object in the
        Segment object. Passes on the remove() call to self.picks.
        """
        self.picks.remove(pick)
        return self

    def select(self, station=None, event=None, channel=None, author=None,
               modes=None, fw=None):
        """
        Return new Segment object only with these picks that match the given
        stats criteria (e.g. all picks with ``station="ADK"``).
        """
        picks = []
        for pick in self:
            # skip trace if any given criterion is not matched
            if station is not None:
                if not fnmatch.fnmatch(pick.stats.station.upper(),
                                       station.upper()):
                    continue
            if event is not None:
                if not fnmatch.fnmatch(pick.stats.event.upper(),
                                       event.upper()):
                    continue
            if channel is not None:
                if not fnmatch.fnmatch(pick.stats.channel.upper(),
                                       channel.upper()):
                    continue
            if author is not None:
                if not fnmatch.fnmatch(pick.stats.author.upper(),
                                       author.upper()):
                    continue
            if modes is not None:
                i = 0
                for mode in pick.stats.modes:
                    if not fnmatch.fnmatch(mode.name.upper(),
                                           modes.upper()):
                        i += 1
                if i == len(pick.stats.modes):
                    continue
            if fw is not None:
                if isinstance(fw, list) and len(fw) == 2:
                    if pick.fw1 < float(fw[0]) or pick.fw2 > float(fw[1]):
                        continue
            picks.append(pick)
        return self.__class__(picks=picks)

    def sort(self, keys=['station', 'fw1', 'fw2',
                         'tw1', 'tw2', 'weight'], reverse=False):
        """
        Sort the picks in the Segment object.

        """
        # check if list
        msg = "keys must be a list of strings. Always available items to " + \
            "sort after: \n'station', 'fw1', 'fw2', 'tw1', 'tw2', 'weight'"
        if not isinstance(keys, list):
            raise TypeError(msg)

        # Loop over all keys in reversed order.
        for _i in keys[::-1]:
            self.picks.sort(key=lambda x: getattr(x, _i), reverse=reverse)

        c = {}
        for i, p in enumerate(self.picks):
            if p.station in ('FURI', 'FUR', 'APE', 'APEZ', 'TRIS', 'TRI'):
                c[p.station] = i

        if 'FURI' in c and 'FUR' in c:
            self.picks = swap(self.picks, c['FURI'], c['FUR'])

        if 'APEZ' in c and 'APE' in c:
            self.picks = swap(self.picks, c['APEZ'], c['APE'])

        if 'TRIS' in c and 'TRI' in c:
            self.picks = swap(self.picks, c['TRIS'], c['TRI'])

        return self

    def write(self, filename, eventfiles=False, overwrite=False, format=None,
              allevents=False, verbose=False):
        """
        filename: will be ignored if eventfiles is set True
        eventfiles: if set True, one segment-file per event is created.
                    The names are based on pick.stats.event
        """
        if len(self.picks) == 0 and overwrite is False:
            return

        # Sort segments by station name, outputs a list of names (segs)
        if eventfiles is True:
            if os.path.isdir(filename):
                segpath = filename
            else:
                segpath = ''
            events = []
            for pick in self:
                if pick.stats.event not in events:
                    events.append(pick.stats.event)
            for ev in events:
                seg = self.select(event=ev)
                seg.sort()
                if verbose:
                    print(ev)
                    print(seg)

                if format == 'pickle':
                    file = os.path.join(segpath, "%s.pickle" % ev)

                elif format == 'segment':
                    fname = "%s.segment%s" % (ev, self[0].stats.channel[-1])
                    file = os.path.join(segpath, fname)

                else:
                    file = os.path.join(segpath, "%s.dat" % ev)

                _write_dat(seg, file)
            return

        if allevents is True:
            with open('allevents', 'w') as fh:
                allevents = []
                for pick in self:
                    allevents.append("%s %s\n" % (pick.stats.channel[-1],
                                     pick.stats.event))
                allevents = list(set(allevents))

                for line in allevents:
                    fh.write(line)
            return
        try:

            # Get filename
            if not overwrite:
                filename = find_unique_name(filename)

            if format is None:
                # try to guess format from file extension
                _, format = os.path.splitext(filename)
                format = format[1:]

            if format == 'pickle':
                write_pickle(self, filename)
            elif format == 'segment':
                fname = "%s.segment%s" % (filename, self[0].stats.channel[-1])
                _write_dat(self, fname)
            elif format in ('db', 'sqlite3', 'sql'):
                _write_db(self, filename, verbose=verbose)
            else:
                _write_dat(self, filename)
        except IOError:
            msg = "\033[91mCan't save file\n"
            msg += "Error message: %s\033[0m" % sys.exc_info()[1]
            print(msg)
        return

    def add_channel(self, channel):
        for p in self:
            p.stats.channel = channel
        return


def _write_dat(segment, filename):
    with open(filename, 'w') as fh:
        # for each stationname 'key' fw and tw is written
        for pick in segment:
            fh.write("%s\n " % str(pick.stats.station))
            fh.write("\t%f\t%f\n" % (pick.fw1, pick.fw2))
            fh.write("\t%f\t%f\n" % (pick.tw1, pick.tw2))
            fh.write("\t%E\n" % pick.weight)


def _write_db(segment, db_path, verbose=False):
    values = {}
    tables = []
    for pick in segment:
        chan = pick.stats.channel
        if chan not in values:
            values[chan] = []
            tables.append(chan)

        if pick.stats.snr is None:
            snr = 'None'
        else:
            snr = pick.stats.snr

        if pick.stats.modes == '':
            modes = ''
        else:
            modes = ','.join(pick.stats.modes.names)

        v = (pick.station, pick.fw1, pick.fw2, pick.tw1, pick.tw2, pick.weight,
             pick.stats.event, modes, snr,
             pick.stats.author, pick.stats.data_origin)

        values[chan].append(v)

    db = sqlite3.connect(db_path)
    c = db.cursor()

    for t in tables:
        if verbose is True:
            print(t)
        schema = "CREATE TABLE if not exists '{table}"
        schema += "' (station text, fw1 real, fw2 real, tw1 real, tw2 real, "
        schema += "weight real, event text, modes text, snr real, "
        schema += "author text, data_origin text)"
        table = schema.format(table=t)

        c.execute(table)
        exe = "INSERT INTO '{table}' VALUES (?,?,?,?,?,?,?,?,?,?,?)"
        exe = exe.format(table=t)
        if verbose is True:
            print(values[t])
        c.executemany(exe, values[t])
    db.commit()
    db.close()
    return
