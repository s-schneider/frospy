from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport
from obspy.core import AttribDict


def init():
    global colors
    # used_station_c = '#e57248'
    # unused_station_c = '#bc8f8f'
    # all_station_c = '#2A98FF'

    colors = AttribDict(
                        ocean_color='#EBEBEB',
                        land_color='#FBFBF2',
                        used_station_c='#e57248',
                        unused_station_c='#bc8f8f',
                        all_station_c='#e57248',
                        marked_station_c='yellow',
                        event_color='#',
                        marker_size=4
                        )
