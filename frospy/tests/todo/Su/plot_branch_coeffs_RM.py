#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:23:08 2019

@author: talavera
"""
import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, basename, splitext
from obspy.core import AttribDict
from collections import OrderedDict
from frospy.postprocessing.uncertainties import uncertainties_calculation
from frospy.postprocessing.AnalyseRuns.plot import misfits_cst, listdir
from frospy.postprocessing import plot
from frospy.postprocessing.plot import inv_summary, cst
from frospy.postprocessing.plot import plot_coeffs_per_mode
from frospy.postprocessing.read import read_inversion_summary
from frospy.postprocessing.misfit import plot_Qf
from frospy.core.splittingfunc import Set
from frospy.core.splittingfunc import loadmodel # python2 branch
#from frospy.splitting.load import loadmodel # python3 master branch
from frospy.core.splittingfunc.read import read_cst, get_cst, read_cst_S20RTS
from frospy.util.read import read_modes_in, get_mode_names
from frospy.util.base import sort_human 
from frospy.core.setup.settings import Setup, get_mode_list
from frospy.util.base import (chunking_list, split_digit_nondigit,
                            max_cc_degrees, max_sc_degrees)
from frospy.core.setup.settings import read as read_setup
from frospy.postprocessing.uncertainties import uncertainties_calculation
from frospy.core.database.write import write_cst
from frospy.core.database.query import get_db_tables, cst_query, db_query
from frospy.core.modes import read as read_modes

##---------------------------------------------------------------------------
master_db_path = '/net/home/talavera/eejit/data/mycst.sqlite3'
#modes_all = read_modes()

sf = Set()         

modes = [
#         ["0S0"],
#         ["1S0-4S2", "1S0","4S2", "0S10"],
         ["2S0-7S2", "2S0","7S2"],
         ["3S0-8S2", "3S0-9S2", "8S2-9S2", "3S0","8S2","9S2"],
         ["4S0-10S2", "4S0-11S2", "10S2-11S2", "4S0","10S2","11S2"],
         ["5S0-13S2", "5S0","13S2"],
         ["6S0"],
         ["7S0"],
         ["8S0"],
         ["9S0"],
         ["11S0-27S2", "11S0", "27S2"],
         ]
db_model = 'GC'
for m in modes:
#    print(m)
    sf += loadmodel(ifile=master_db_path, modes=m,
                   name=db_model, damp='0', db_model=db_model)    

modes = [
#         "0S0","1S0",
         "2S0","3S0","4S0","5S0","6S0","7S0","8S0","9S0",
#         "11S0","27S2",
        ]
db_model = 'SC'
for m in modes:
#    mode = modes_all.select(name=m)[0]
#    print(mode.name, mode.freq*1000, mode.Q)

    #    print(m)
    sf += loadmodel(ifile=master_db_path, modes=m,
                   name=db_model, damp='0', db_model=db_model)



## center frequency
#SF = plot_coeffs_per_mode(SF_in=sf, 
#                          label1 = "GC", 
#                          label2='SC', 
#                          mode='sc', 
#                          model=['$\phi$=1.20','$\phi$=1.10','$\phi$=1.04',
#                                 '$\phi$=0.96','$\phi$=0.90','$\phi$=0.80',
#                                 "S20RTS",], #
#                          plot_f_center=True, 
#                          degree = 0, 
#                          ordering=["l", "sens"],
#                          l=0,
#                          cmap='GreensBluesBlack', 
#                          rot=0,
#                          markersize = 6,
#                          color1="r",
#                          color2="k",
#                          fig_size = (7,4),
##                          ylim=[-3.5,9],
##                          savefig=True,
##                          verbose=True,
#                          )
#
#
## center frequency
#SF = plot_coeffs_per_mode(SF_in=sf, 
#                          label1 = "GC", 
#                          label2='SC', 
#                          mode='sc', 
#                          model=["GD","QM1","DE","HT", "REM","Sumatra"], #["AD","GLW","REM"],#
##                          model=["QM1","DE","HT", "REM","Sumatra"], #["AD","GLW","REM"],#
#                          plot_f_center=True, 
#                          degree = 0, 
#                          ordering=["l", "sens"],
#                          l=0,
#                          cmap='grey', 
#                          rot=0,
#                          spacing = True,
#                          markersize = 6,
#                          color1="r",
#                          color2="k",
#                          fig_abc="a",
#                          fig_size = (8,3.5),
#                          yaxs_split=[-25,-8,-3.5,10,],
##                          ylim=[-3.5,9],
##                          savefig=True,
##                          verbose=True,
##                          legend_show=False,
#                          )
#
## quality factor
#SF = plot_coeffs_per_mode(SF_in=sf,  
#                          label1 = "GC", 
#                          label2='SC', 
#                          mode='sc', 
#                          model=["QM1"],
##                          model=["GD","QM1","DE","HT", "REM","Sumatra"], #["AD","GLW","REM"],#
##                          model=["QM1","DE","HT", "REM","Sumatra"], #["AD","GLW","REM"],#
#                          plot_Q=True, kind = "dst", 
#                          degree = 0, 
#                          l=0,
#                          ordering=["l", "sens"],
#                          cmap='grey',
#                          rot=0,
#                          spacing = True,
#                          markersize = 6,
#                          color1="r",
#                          color2="k",
#                          fig_abc="b",
#                          fig_size = (8.1,3.5),
##                          ylim=[-250,1100],
##                          ylim=[-50,1000],
##                          savefig=True,
##                          verbose=True,
##                          legend_show=True,
#                          )

# degree two # if with location of GLW to AD
SF = plot_coeffs_per_mode(SF_in=sf,  label1 = "GC",
                          label2='SC', 
                          mode='sc', 
                          model=["Wh","Ro","Tr","BT"],#
#                          model=["S20RTS","AD","GLW",],#
#                          model=["S20RTS","Tr","AD","GLW",],#
                          degree = 2, 
                          order=0,
                          l=2,  
                          ordering=["sens", "l"], 
                          cmap='grey',
                          rot=0, 
                          spacing = True,
                          markersize = 6.,
                          fig_abc="a",
#                          yaxs_split=[-50,140,186,440], #all IC models
#                          yaxs_split=[-18,60,150,500], #Tr
#                          yaxs_split=[-18,80,100,150], #Wh
#                          yaxs_split=[-18,40,60,70], #Ro
#                          yaxs_split=[-18,65,100,150], #BT
#                          ylim=[-10,40],
#                          verbose=True,
#                          savefig=True,
                          )

## degree four
SF = plot_coeffs_per_mode(SF_in=sf,   label1 ='GC', 
                          label2='SC', 
                          mode='sc', 
                          model=["Wh","Ro","Tr","BT"],#
#                          model=["S20RTS","AD","GLW",],#, #'S20RTS+CRUST+BT'],  
#                          model=["S20RTS","Tr","AD","GLW",],#, #'S20RTS+CRUST+BT'],  
                          degree = 4, 
                          order=0,
                          l=2,  
                          ordering=["sens", "l"], 
                          cmap='grey',
                          rot=0, 
                          spacing = True,
                          markersize = 6.,
                          fig_abc="b",
#                          yaxs_split=[-10,25,50,140],#Tr
#                          yaxs_split=None,#Wh
#                          yaxs_split=[-160,-100,-18,40],#Ro
#                          ylim=[-10,24],
#                          legend_show=False,
#                          savefig=True,
                          )

## CC degrees
SF = plot_coeffs_per_mode(SF_in=sf,   label1 ='GC', 
                          mode='cc', 
#                          model=['S20RTS'], 
                          model=["Wh","Ro","Tr","BT"],#
#                          model=["S20RTS"], #'S20RTS+CRUST+BT'],  
#                          model=["S20RTS","Tr",], #'S20RTS+CRUST+BT'],  
                          degree=2, 
                          order=0,
                          l=[0,2],
                          ordering=["l","sens"], 
                          cmap='grey', 
                          spacing = True,
                          markersize = 6.,
                          rot=0,
#                          verbose=True,
                          fig_abc="c",
#                          yaxs_split=[-60,10,35,38],#Ro
#                          legend_show=False,
#                          fig_size = (8,5),
#                          savefig=True,
                          )
