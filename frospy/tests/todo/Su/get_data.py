h#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 16:36:40 2018

@author: talavera
"""

import numpy as np
import matplotlib.pyplot as plt
from frospy.util.read import read_st

data_raw='/data/talavera/data/060994A.13dy.ahx'
st_raw = read_st(data_raw, 'ah')
st_raw.sort(['station'])

data='/data/talavera/data/060994A.13dy.noglitch.ahx'
st_tmp = read_st(data, 'ah')
st_tmp.sort(['station'])

data_new='/data/talavera/data/060994A.13dy.new.noglitch.ahx'
st_tmp.write(data_new, format='ah')

#from frospy.util.read import read_st
#event='060994A'
#data='300hrs/raw/VHZ/%s/%s.v3.ahx'%(event,event)
#
#st_tmp = read_st(data, 'ah')
#st_tmp
#st_tmp.sort(['station'])
#st_tmp.write(data, format='ah')
