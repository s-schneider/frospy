#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 18:26:24 2019

@author: talavera
"""

import os
import sys
import glob
import numpy as np
from os.path import join, basename, splitext
from frospy.core.segment import read
from obspy.core import AttribDict
from frospy.core.database.query import db_query
from frospy.core.database.query import get_db_tables
from frospy.core.segment import read

cdir = '/net/home/talavera/eejit/splitting/'
mode  = '03s00-08s02-09s02'
event ='%sevents'%'good'
cpl   ='cst+d00+cc' #c00+d00 cst+d00
mzero ='08s02-09s02' #romanowicz cst+d00_best_result

db_path = join(cdir,mode,event,cpl,mzero, "inversion_out", "inversion.db")
damp=0.05
iteration=15
seg_suffix = 'dat',
#sets_name = "08s02-09s02"
sets_name = "03s00"

allevents = set(db_query(db_path, sets_name,
                condition_dict={'iter': iteration, 'damp': damp},
                select='event'))


allmisfits = []
for e in list(allevents):
    e=map(str, e)[0]
    if e !='122604Q':
        misfits = set(db_query(db_path, sets_name,
                               condition_dict={'iter': iteration, 'damp': damp,
                                               'event': e},
                                               select='misfit'))
        allmisfits.extend(list(misfits))
        print e, len(misfits), np.mean(list(misfits))
seg = len(allmisfits)
misfit = np.mean(list(allmisfits))
print(sets_name, seg, misfit)
#
#('01s00', 241, 0.09825035684821577) 0.1029	241
#('03s00', 332, 0.39368779346536137) 0.4230	349
#('04s00', 487, 0.5001769146484395)  0.1685	302
#('05s00', 146, 0.2131315034609589) 0.2252	146



import os
import sys
import glob
import numpy as np
from os.path import join, basename, splitext
from frospy.core.segment import read
from obspy.core import AttribDict
from frospy.core.database.query import db_query
from frospy.core.database.query import get_db_tables
from frospy.core.segment import read

cdir = '/net/home/talavera/eejit/splitting/'
mode  = '00t09'
event ='cst+dst'#'%sevents'%'all'
cpl   ='prem_TR' #c00+d00 cst+d00
mzero ='cst18+dst12' #romanowicz cst+d00_best_result
#mode  = '03s00'
#event ='%sevents'%'all'
#cpl   ='self' #c00+d00 cst+d00
#mzero ='prem' #romanowicz cst+d00_best_result


db_path = join(cdir,mode,event,cpl,mzero, "inversion_out", "inversion.db")
damp=0.01
iteration=10
sets_name = "VHT"

allevents = set(db_query(db_path, sets_name,
                condition_dict={'iter': iteration, 'damp': damp},
                select='event'))

print("tot events", len(allevents))
allmisfits = []
for e in list(allevents):
    e=map(str, e)[0]
#    if e !='122604A':
    misfits = set(db_query(db_path, sets_name,
                           condition_dict={'iter': iteration, 'damp': damp,
                                           'event': e},
                                           select='misfit'))
    allmisfits.extend(list(misfits))
    print(e,len(misfits),np.mean(list(misfits)))

seg = len(allmisfits)
misfit = np.mean(list(allmisfits))
print(sets_name, seg, misfit)
