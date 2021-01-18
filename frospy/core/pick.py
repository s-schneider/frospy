# -*- coding: utf-8 -*-
"""
Module for handling nmPy pick objects.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from copy import deepcopy
from obspy.core.trace import Trace
from obspy.core.util import AttribDict
from frospy.core.modes import Modes, Mode
from getpass import getuser


class Stats(AttribDict):
    defaults = {
        'modes': '',
        'station': '',
        'event': '',
        'channel': '',
        'snr': None,
        'author': '',
        'data_origin': '',
        'misfit': None
    }

    _refresh_keys = {'modes', 'station', 'event', 'channel', 'snr', 'author',
                     'data_origin', 'misfit'}

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
            if key in ['station', 'event', 'channel', 'author',
                       'data_origin']:
                value = str(value)
            elif key == 'modes':
                if isinstance(value, Modes) or isinstance(value, Mode):
                    value = value
                elif type(value) is list:
                    value = [str(val) for val in value]
                else:
                    value = str(value)

            # equivalent to AttribDict.__setitem__(self, key, value)
            super(Stats, self).__setitem__(key, value)

    __setattr__ = __setitem__

    def __str__(self, extended=False):
        """
        Return better readable string representation of Stats object
        """
        _pretty_str = '%s\t|' % (self.station)
        _pretty_str += ' channel: %s | event: %s |' % (self.channel,
                                                       self.event)
        if isinstance(self.modes, Modes):
            mode_list = []
            for mod in self.modes:
                mode_list.append(mod.name)
            _pretty_str += ' mode(s): %s' % (mode_list)
        else:
            _pretty_str += ' mode(s): %s' % (self.modes)
        if extended is True:
            _pretty_str += ' | author: %s |' % (self.author)
            _pretty_str += ' data origin: %s' % (self.data_origin)
        return _pretty_str

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))


class Pick(object):
    """
    An object containing values of a frequency window, time window and an
    assigned weight corresponding to the amplitude spectrum of a timeseries.

    :type station: string
    :param station: Code of the recording station
    :type fw1, fw2: float
    :param fw1, fw2: Corner points of the given frequency interval in mHz
    :type tw1, tw2: float
    :param tw1, tw2: Corner points of the given time interval in hrs
    :type weight: float
    :param weight: Weight of this pick for later use in inversion
    """
    def __init__(self, header=None, fw1=0, fw2=0, tw1=0, tw2=0, weight=0):
        # property of pick:
        if header is None:
            header = {}
        elif type(header) is Trace:
            header_tmp = {'station': header.stats.station,
                          'channel': header.stats.channel}
            header = header_tmp

        self.stats = Stats(header)
        self.stats.author = getuser()
        self.station = self.stats.station
        self.fw1 = float(fw1)
        self.fw2 = float(fw2)
        self.tw1 = float(tw1)
        self.tw2 = float(tw2)
        self.weight = float(weight)

    def __str__(self):
        """
        Return better readable string representation of Pick object.
        """
        _pretty_str = '%s\t|' % (self.stats.station)
        _pretty_str += ' fw1: %.3f fw2: %.3f |' % (self.fw1, self.fw2)
        _pretty_str += ' tw1: %.3f tw2: %.3f |' % (self.tw1, self.tw2)
        _pretty_str += ' weight: %.2e' % (self.weight)
        return _pretty_str

    def __eq__(self, other):
        return self.stats.station == other.stats.station and \
               self.fw1 == other.fw1 and self.fw2 == other.fw2 and \
               self.tw1 == other.tw1 and self.tw2 == other.tw2 and \
               self.weight == other.weight

    def __hash__(self):
        return (hash(('station', self.stats.station,
                      'fw1', self.fw1,
                      'fw2', self.fw2,
                      'tw1', self.tw1,
                      'tw2', self.tw2,
                      'weight', self.weight)))

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def copy(self):
        return deepcopy(self)
