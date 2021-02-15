# -*- coding: utf-8 -*-
"""
Converts pickled segment files to .dat files compatable with synseis.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from frospy.core.segments import read
import sys

if len(sys.argv) < 2:
    print("USAGE python segment_pickle_to_dat.py pickle_file")
    sys.exit()

filename = sys.argv[1]
segment = read(filename)

cmt_id = segment[0].stats.event
segment.write("%s.dat", cmt_id)
