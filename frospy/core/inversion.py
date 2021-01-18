# -*- coding: utf-8 -*-
from obspy.core import AttribDict

"""
Module for handling nmPy cst_inv objects.
Potential use of ASDF for file storage

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""


class Stats(AttribDict):
    """
    Stats object for Spectrum class
    """
    # defaults = {
    #             'allmisfits': []
    #             }
    #
    # _refresh_keys = {'allmisfits'}

    def __init__(self, header={}):
        """
        """
        # super() is used to initialize AttribDict within Stats(),
        # which sets the defaults as initial dictionary
        # equivalent to AttribDict.__init__(self, header)
        super(Stats, self).__init__(header)


class Inversion(object):
    """
    This class contains all relevant data and variables created for and
    resulting from a cst or fsi inversion, listed as follows:

    Spectrum object:
        catalog of used events :obspy.core.Catalog:
        streams of seismic data (per event) :obspy.core.Stream:
        streams of synthtic data (per event) :obspy.core.Stream:
        inventory of used stations (per event) :obspy.core.Inventory:
        Segments (per event) :frospy.core.Segment:
        Modes that were inverted for :frospy.core.modes.Modes:

    SplittinFunc or FSI object:
        the result of the cst inversion

    Stats object:
        Used values, used as input for the .conf files:
        damping, iterations, misfit, setup of startmodel etc
    """
    def __init__(self, stats={}):
        self.stats = Stats(stats)
        return

    def append_stats(self, stats):
        if type(stats) in [dict, AttribDict]:
            for cmt, values in stats.items():

                if cmt in self.stats:
                    self.stats[cmt]['n'] = values['n']
                else:
                    self.stats[cmt] = {'misfit': {}, 'n': values['n']}
                for it, mf in values['misfit'].items():
                    self.stats[cmt].misfit[it] = mf
        return

    def append_allmisfits(self, misfits):
        self.allmisfits = misfits
        return
