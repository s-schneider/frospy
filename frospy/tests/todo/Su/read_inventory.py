#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 09:55:46 2020

@author: sujania
"""

from frospy.util.data_request import process_downloader_data
from frospy.util.read import read_inventory as read_inv
from frospy.util.read import read_std_cat

inv1 = read_inv('/Users/sujania/code/nmPy/nmpy/data/AD/arwens_stations.xml')
inv2 = read_inv('/Users/sujania/code/nmPy/nmpy/data/AD/stations.xml')

for i in inv1:
    print(("%s %03i %s")%(i.code, i.total_number_of_stations, i.description))
    
for i in inv2:
    print(("%s %03i %s")%(i.code, i.total_number_of_stations, i.description))