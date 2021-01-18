#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 10:39:00 2018

@author: talavera
"""

from frospy.preprocessing.segment_picker import autopick_event
import os
import glob


def pick_segments(run_dir, data_path, events, comps, modes=None,
                  fw=None, tw=None, max_mf=None,
                  verbose=False, overwrite=False):
    """
    from frospy.preprocessing.pick_segments import pick_segments
    from frospy.core.modes import read as read_modes
    events = ['060994A', '100494B', '031111B', '022710A']
    data_path = '//nfs/stig/simons/alldata'
    comps = ['VHT']
    min_snr=2
    fw = [1.715, 1.75]
    tw = [3,60]
    modes = ['00t11-00s10']
    run_dir = '//nfs/stig/simons/splitting/cross_coupling'
    overwrite = False
    from frospy.preprocessing.pick_segments import pick_segments

    pick_segments(run_dir, data_path, events, comps, modes,
                  fw, tw, min_snr, verbose=True, overwrite=overwrite)

    """
    N = len(events)
    for i, event in enumerate(events):
        msg = '\n\x1b[6;30;42m\n\nCurrent event [%i/%i]: %s\n\x1b[0m'
        print(msg % (i+1, N, event))
        autopick_event(data_path, run_dir, comps, modes,
                       min_snr, fw=fw, event=event, tw=tw,
                       verbose=verbose, overwrite=overwrite)

    return

events = ['011307A', '022710A', '030994E', '031111B', 
          '032598B', '032805D', '050306F', '052413A', 
          '053015A', '060994A', '061796A', '070508A', 
          '073095A', '100494B', '111506F', '122604A']
data_path = '//nfs/stig/simons/alldata'
comps = ['VHZ']
min_snr=1
fw = [7.400000, 7.438000]
tw = [10,55]
modes = ['8s0']
overwrite = True
run_dir = '//nfs/stig/talavera/radial-inversion/08s00/segments/08s00-all'

pick_segments(run_dir, data_path, events, comps, modes,
              fw, tw, min_snr, verbose=True, overwrite=overwrite)