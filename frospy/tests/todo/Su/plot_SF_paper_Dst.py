#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:06:26 2019

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
from frospy.postprocessing.plot import plot_coeffs_per_mode

master_db_path = '/net/home/talavera/eejit/data/mycst.sqlite3'  

modes = [
#         ["0S2"],
#         ["0S5"],
#         ["0S6"],
#         ["0S7"],
#         ["1S4"],
#         ["1S7"],
#         ["1S8"],
#         ["1S9"],
#         ["1S10"],
#         ["2S6"],
         ["2S12"],
#         ["2S13"],
#         ["3S6"],
#         ["3S9"],
#         ["0T6"],
#         ["0T8"],
#         ["0T9"],
#         ["1T2"],
#         ["1T5"],
#         ["2T3"],
#         ["2T7"],
         ]

scut = [
#         [2,2],
#         [8,6],
#         [10,4],
#         [12,6],
#         [6,4],
#         [10,6],
#         [12,6],
#         [12,6],
#         [12,4],
#         [10,6],
         [12,10],
#         [12,2],
#         [6,2],
#         [12,4],
#         [6,2],
#         [8,2],
#         [8,4],
#         [2,2],
#         [8,2],
#         [6,4],
#         [6,2],
         ]

db_model = 'dst'
model = "S20RTS"
savefigure=False
for m,smax in zip(modes,scut):
#    print(m)
    sf = Set()
    sf += loadmodel(ifile=master_db_path, modes=m,
                   name=db_model, damp='0', db_model=db_model)  
    SF = plot_coeffs_per_mode(SF_in=sf,
                          degree=2,
                          kind="dst",
                          model=[model],
                          )
    SF[0].stats.name = "Our Inversion"
    SF[1].stats.name = model
    if m==["0S6"] or m==["2S12"]:
        p = SF.plot_map(
                    smax=smax[0],
                    R=-0.2,
                    kind="cst",
                    kind_in_title=False, 
        #            fig_size=[7,1.5],
                    show_title=False,
                    legend_show=True, 
    #                fig_abc='b',
#                    vmax=int(6),
#                    vmin=int(-11),
                    savefigure=savefigure,
                    )
    else:
        p = SF.plot_map(
                    R=-0.2,
                    smax=smax[0],
                    kind_in_title=False, 
        #            fig_size=[7,1.5],
                    show_title=False,
                    legend_show=False, 
    #                fig_abc='a',
        #            vmax=3.5,
        #            vmin=-4,
                    savefigure=savefigure,
                    )
        
    if m==["1T2"]:
        p = SF.plot_map(
                    smax=smax[1],
                    R=-0.2,
                    kind="dst",
                    kind_in_title=False, 
        #            fig_size=[7,1.5],
                    show_title=False,
                    legend_show=False, 
    #                fig_abc='b',
                    vmax=0.1,
                    vmin=-0.1,
                    savefigure=savefigure,
                    )
    elif m==["0T6"]:
        p = SF.plot_map(
                    smax=smax[1],
                    R=-0.2,
                    kind="dst",
                    kind_in_title=False, 
        #            fig_size=[7,1.5],
                    show_title=False,
                    legend_show=False, 
    #                fig_abc='b',
                    vmax=0.3,
                    vmin=-0.3,
                    savefigure=savefigure,
                    )
    else:
            p = SF.plot_map(
                smax=smax[1],
                R=-0.2,
                kind="dst",
                kind_in_title=False, 
    #            fig_size=[7,1.5],
                show_title=False,
                legend_show=False, 
#                fig_abc='b',
    #            vmax=3.5,
    #            vmin=-4,
                savefigure=savefigure,
                )