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
#sf = Set()   
# 
#dhome = "/net/home/talavera/eejit/splitting/"
#mdir = [ # w/ 122604A & 052389A
#"00s02/cst/prem/cst4/inversion_out/mnew-it10-d1.00e-03.dat",
#"00s05/cst/prem/cst10/inversion_out_4/mnew-it10-d1.00e-02.dat",
#"00s06/cst/prem/cst12/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"00s07/cst/prem/cst14/inversion_out_4/mnew-it10-d1.00e-02.dat",
#"01s04/cst/prem/cst8/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"01s07/cst/prem/cst14/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"01s08/cst/prem/cst16/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"01s09/cst/prem/cst18/inversion_out_6/mnew-it10-d1.00e-02.dat",
#"01s10/cst/prem/cst20/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"02s06/cst/prem/cst12/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"02s12/cst/prem/cst20/inversion_out_3/mnew-it10-d1.00e-02.dat",
#"02s13/cst/prem/cst20/inversion_out_4/mnew-it10-d1.00e-02.dat",
#"03s06/cst/prem/cst12/inversion_out_4/mnew-it10-d1.00e-02.dat",
#"03s09/cst/prem/cst18/inversion_out_3/mnew-it10-d1.00e-02.dat",
##"04s09/cst/prem/cst18/inversion_out_4/mnew-it10-d1.00e-02.dat",
##"05s03/cst/prem/cst6/inversion_out_4/mnew-it10-d1.00e-02.dat",
##"05s07/cst/prem/cst14/inversion_out_3/mnew-it10-d1.00e-02.dat",
##"05s08/cst/prem/cst16/inversion_out_5/mnew-it10-d1.00e-02.dat",
##"06s10/cst/s20rts/cst20/inversion_out_3/mnew-it10-d1.00e-02.dat",
##"07s05/cst/prem/cst10/inversion_out_6/mnew-it10-d1.00e-02.dat",
##"07s07/cst/s20rts/cst14/inversion_out_5/mnew-it10-d1.00e-02.dat",
##"08s06/cst/s20rts/cst12/inversion_out_6/mnew-it10-d1.00e-03.dat",
##"09s11/cst/s20rts/cst20/inversion_out_3/mnew-it10-d1.00e-03.dat",
##"11s09/cst/s20rts/cst18/inversion_out_3/mnew-it10-d1.00e-03.dat",
##"12s06/cst/s20rts_AD/cst12/inversion_out_2/mnew-it10-d1.00e-03.dat",
##"12s07/cst/s20rts_AD/cst14/inversion_out_2/mnew-it10-d1.00e-04.dat",
##"12s13/cst/s20rts_AD/cst20/inversion_out_2/mnew-it10-d1.00e-03.dat",
#]

#mdir = [ # w/o 122604A & 052389A
#"00s02/cst/prem/cst4/inversion_out/mnew-it10-d1.00e-03.dat",
#"00s05/cst/prem/cst10/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s06/cst/prem/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s07/cst/prem/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s04/cst/prem/cst8/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s07/cst/prem/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s08/cst/prem/cst16/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s09/cst/prem/cst18/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s10/cst/prem/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s06/cst/prem/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s12/cst/prem/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s13/cst/prem/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s06/cst/prem/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s09/cst/prem/cst18/inversion_out/mnew-it10-d1.00e-02.dat",
#"04s09/cst/prem/cst18/inversion_out/mnew-it10-d1.00e-02.dat",
#"05s03/cst/prem/cst6/inversion_out/mnew-it10-d1.00e-02.dat",
#"05s07/cst/prem/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"05s08/cst/prem/cst16/inversion_out/mnew-it10-d1.00e-02.dat",
#"06s10/cst/s20rts/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"07s05/cst/prem/cst10/inversion_out/mnew-it10-d1.00e-02.dat",
#"07s07/cst/s20rts/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"08s06/cst/s20rts/cst12/inversion_out/mnew-it10-d1.00e-03.dat",
#"09s11/cst/s20rts/cst20/inversion_out/mnew-it10-d1.00e-03.dat",
#"11s09/cst/s20rts/cst18/inversion_out/mnew-it10-d1.00e-03.dat",
#"12s06/cst/s20rts_AD/cst12/inversion_out/mnew-it10-d1.00e-03.dat",
#"12s07/cst/s20rts_AD/cst14/inversion_out/mnew-it10-d1.00e-04.dat",
#"12s13/cst/s20rts_AD/cst20/inversion_out/mnew-it10-d1.00e-03.dat",
#]

#mdir = [ # AD
#"00s02/cst/prem_AD/cst4/inversion_out/mnew-it10-d1.00e-03.dat",
#"00s05/cst/prem_AD/cst10/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s06/cst/prem_AD/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s07/cst/prem_AD/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s04/cst/prem_AD/cst8/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s07/cst/prem_AD/cst14/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s08/cst/prem_AD/cst16/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s09/cst/prem_AD/cst18/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s10/cst/prem_AD/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s06/cst/prem_AD/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s12/cst/prem_AD/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s13/cst/prem_AD/cst20/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s06/cst/prem_AD/cst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s09/cst/prem_AD/cst18/inversion_out/mnew-it10-d1.00e-02.dat",
#]

#ifiles = [glob.glob(join(dhome,m))[0] for m in mdir]
#for f in ifiles:
#    r = f.split('/')[0:-2]
#    main = join('/',*r)
#    setup = join(main,'setup.json')
#    print(setup)
#    setup = read_setup(setup)
#    damp  = setup.damping
#    sf += loadmodel(setup=setup, ifile = f, 
#                    name ='cst', damp = damp)
    
    
#mdir = [
#"00s02/cst+dst/prem_Z/cst4+dst2/inversion_out/mnew-it10-d1.00e-03.dat",
#"00s05/cst+dst/prem_TRZ/cst10+dst10/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s06/cst+dst/prem_TRZ/cst12+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"00s07/cst+dst/prem_TRZ/cst14+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s04/cst+dst/prem_TRZ/cst8+dst8/inversion_out/mnew-it10-d1.00e-01.dat",
#"01s07/cst+dst/prem_TZ/cst14+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s08/cst+dst/prem_TRZ/cst16+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s09/cst+dst/prem_TRZ/cst18+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"01s10/cst+dst/prem_TZ/cst20+dst12/inversion_out/mnew-it10-d1.00e-03.dat",
#"02s06/cst+dst/prem_TRZ/cst12+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s12/cst+dst/prem_TRZ/cst20+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s13/cst+dst/prem_TRZ/cst20+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s06/cst+dst/prem_TRZ/cst12+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
#"03s09/cst+dst/prem_TRZ/cst18+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
##"04s09/cst+dst/prem_TRZ/cst18+dst12/inversion_out/mnew-it10-d1.00e-03.dat",
##"05s03/cst+dst/prem_TRZ/cst6+dst6/inversion_out/mnew-it10-d1.00e-02.dat",
##"05s07/cst+dst/prem_TRZ/cst14+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
##"05s08/cst+dst/prem_TRZ/cst16+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
##"06s10/cst+dst/prem_TRZ/cst20+dst12/inversion_out/mnew-it10-d1.00e-03.dat",
##"07s05/cst+dst/prem_TRZ/cst10+dst10/inversion_out/mnew-it10-d1.00e-02.dat",
##"07s07/cst+dst/prem_TRZ/cst14+dst12/inversion_out/mnew-it10-d1.00e-02.dat",
##"08s06/cst+dst/prem_TRZ/cst12+dst12/inversion_out/mnew-it10-d1.00e-03.dat",
##"09s11/cst+dst/prem_TRZ/cst20+dst12/inversion_out/mnew-it10-d1.00e-03.dat",
##"11s09/cst+dst/prem_Z/cst18+dst2/inversion_out/mnew-it10-d1.00e-03.dat",
##"12s06/cst+dst/prem_AD/cst12+dst6/inversion_out/mnew-it10-d1.00e-03.dat",
##"12s07/cst+dst/prem_AD/cst14+dst6/inversion_out/mnew-it10-d1.00e-03.dat",
##"12s13/cst+dst/prem_AD/cst20+dst6/inversion_out/mnew-it10-d1.00e-04.dat",
#]
#
#ifiles = [glob.glob(join(dhome,m))[0] for m in mdir]
#for f in ifiles:
#    r = f.split('/')[0:-2]
#    main = join('/',*r)
#    setup = join(main,'setup.json')
#    print(setup)
#    setup = read_setup(setup)
#    damp  = setup.damping
#    sf += loadmodel(setup=setup, ifile = f, 
#                    name ='cst+dst', damp = damp)


master_db_path = '/net/home/talavera/eejit/data/mycst.sqlite3'
#modes_all = read_modes()

sf = Set()         

modes = [
         ["2S6"], #UM
         ["2S12"], #UM
#         ["2S13"], #UM
         ["3S6"], #UM
         ["3S9"], #UM
         ["0T6"], #UM
#         ["0T8"], #UM
#         ["0T9"], #UM
         ["1T2"], #UM
#         ["0S2"], #LM
#         ["0S5"], #LM
#         ["0S6"], #LM
#         ["0S7"], #LM
#         ["1S4"], #LM
#         ["1S7"], #LM
#         ["1S8"], #LM
#         ["1S9"], #LM
#         ["1S10"], #LM
#         ["1T5"], #LM
#         ["2T3"], #LM
#         ["2T7"], #LM
         ]
db_model = 'dst'
for m in modes:
#    print(m)
    sf += loadmodel(ifile=master_db_path, modes=m,
                   name=db_model, damp='0', db_model=db_model)    

# Cst's

# degree one
SF = plot_coeffs_per_mode(SF_in=sf,  
#                          label1 = 'cst',
#                          label2 = 'cst+dst',
                          label1 = 'cst+dst',
#                          color1="k",
#                          color2="r",
                          color1="r",
                          mode='sc', 
                          kind="cst",
                          model=["AD",'S20RTS'],
                          degree = 0, 
#                          order=0,
#                          l=2,  
#                          n=2,
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(16,4),
#                          fig_abc="a",
#                          ylim=[-10,40],
#                          verbose=True,
#                          savefig=True,
                          )

# degree two
SF = plot_coeffs_per_mode(SF_in=sf,  
#                          label1 = 'cst',
#                          label2 = 'cst+dst',
                          label1 = 'cst+dst',
#                          color1="k",
#                          color2="r",
                          color1="r",
                          mode='sc', 
                          kind="cst",
                          model=["AD",'S20RTS'],
                          degree = 2, 
                          R=-0.2,
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,8),
#                          fig_abc="a",
#                          ylim=[-10,40],
#                          verbose=True,
                          savefig=True,
                          )

# degree four
SF = plot_coeffs_per_mode(SF_in=sf,   
#                          label1 = 'cst',
#                          label2 = 'cst+dst',
                          label1 = 'cst+dst',
#                          color1="k",
#                          color2="r",
                          color1="r",
                          mode='sc', 
                          kind="cst",
                          model=["AD",'S20RTS'], #'S20RTS+CRUST+BT'],  
                          degree = 4, 
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,12),
#                          fig_abc="b",
#                          ylim=[-10,20],
#                          legend_show=False,
                          savefig=True,
                          )


# degree six
SF = plot_coeffs_per_mode(SF_in=sf,   
#                          label1 = 'cst',
#                          label2 = 'cst+dst',
                          label1 = 'cst+dst',
#                          color1="k",
#                          color2="r",
                          color1="r",
                          mode='sc', 
                          kind="cst",
                          model=["AD",'S20RTS'], #'S20RTS+CRUST+BT'],  
                          degree = 6, 
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,16),
#                          fig_abc="b",
#                          ylim=[-10,20],
#                          legend_show=False,
                          savefig=True,
                          )

# Dst's
#degree zero
SF = plot_coeffs_per_mode(SF_in=sf,  
#                          label1 = 'cst',
#                          label2 = 'cst+dst',
                          label1 = 'cst+dst',
#                          color1="k",
#                          color2="r",
                          color1="r",
                          mode='sc', 
                          kind="dst",
#                          plot_Q=True, 
                          model=["AD"],
                          degree = 0, 
#                          order=0,
#                          l=2,  
#                          n=2,
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(16,4),
#                          fig_abc="a",
#                          ylim=[-10,40],
#                          verbose=True,
#                          savefig=True,
                          )
# degree two
SF = plot_coeffs_per_mode(SF_in=sf,  
#                          label1 = 'cst+dst',
                          label1 = 'Our Inversion',
                          color1="r",
                          mode='sc', 
                          kind="dst",
                          model=["S20RTS"],
                          degree = 2, 
                          R=-0.2,
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,8),
#                          fig_abc="a",
#                          ylim=[-1.2,1],
#                          verbose=True,
#                          savefig=True,
                          )

# degree four
SF = plot_coeffs_per_mode(SF_in=sf,   label1 = 'cst+dst',
                          color1="r",
                          mode='sc', 
                          kind="dst",
                          model=["S20RTS"], #'S20RTS+CRUST+BT'],  
                          degree = 4, 
                          R=-0.2,
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,12),
#                          fig_abc="b",
#                          ylim=[-10,20],
#                          legend_show=False,
#                          savefig=True,
                          )


# degree six
SF = plot_coeffs_per_mode(SF_in=sf,   label1 = 'cst+dst',
                          color1="r",
                          mode='sc', 
                          kind="dst",
                          model=["S20RTS"], #'S20RTS+CRUST+BT'],  
                          degree = 6, 
                          R=-0.2,
#                          order=0,
#                          l=2,  
#                          ordering=["freq"], 
                          ordering=["n", "l"], 
                          cmap='grey',
#                          rot=0, 
#                          spacing = True,
                          markersize = 6.,
                          fig_size=(12,16),
#                          fig_abc="b",
#                          ylim=[-10,20],
#                          legend_show=False,
                          savefig=True,
                          )

#SF = plot_coeffs_per_mode(SF=sf,  label1 = "CC",
#                          mode='cc', model=["AD", "GLW", 'S20RTS+CRUST+BT'], 
#                          degree=4, ordering="l")

#dhome = "/net/home/talavera/eejit/splitting/"
#mdir = [
#"01s00-04s02-00s10/*/cst+d00+cc/prem/inversion_out/mnew-it10-d1.00e-02.dat",
#"02s00-*/*/cst+d00+cc/cst_02s00_c20=15/inversion_out/mnew-it10-d1.00e-03.dat",
#"08s02-09s02/*/cst+d00+cc/begtramp/inversion_out_5/mnew-it20-d1.00e-02.dat",
#"04s00-*/*/cst+d00+cc/10s02-11s02/inversion_out_2/mnew-it10-d1.00e-02.dat",
#"05s00-*/*/cst+d00+cc/cst_05s00_13s02/inversion_out/mnew-it10-d1.00e-04.dat",
#"11s00-*/*/cst+d00+cc/cst_27s02+c00=20/inversion_out/mnew-it20-d1.00e-03.dat",
#]
#
#ifiles = [glob.glob(join(dhome,m))[0] for m in mdir]
#
#mdir = [
#"01s00-04s02-00s10/*/cst+d00/prem/inversion_out/mnew-it10-d1.00e-03.dat",
#"02s00-*/allevents/cst+d00/cst_02s00/inversion_out/mnew-it10-d1.00e-03.dat",
#"03s00-*/*/cst+d00/08s02-09s02/inversion_out/mnew-it10-d1.00e-02.dat",
#"04s00-*/*/cst+d00/10s02-11s02/inversion_out/mnew-it10-d1.00e-02.dat",
#"05s00-*/*/cst+d00/cst_05s00_13s02/inversion_out/mnew-it10-d1.00e-04.dat",
#]
#
#ifiles2 = [glob.glob(join(dhome,m))[0] for m in mdir]
#
## center frequency
#SF = plot_coeffs_per_mode(ifiles, label1 = "CC",
#                          ifiles2=ifiles2, label2="ellip CC", 
#                          mode='sc', model=["AD","HT", "QM1", "REM"], 
#                          plot_f_center=True, 
#                          degree = 0, ordering=["sens", "l"], l=2,
#                          cmap='grey', rot=0,
##                          savefig=True,
#                          )
#
## quality factor
#SF = plot_coeffs_per_mode(ifiles, label1 = "CC", kind = "dst",
#                          ifiles2=ifiles2, label2="ellip CC", 
#                          mode='sc', model=["AD","HT", "QM1", "REM", "DE"], 
#                          plot_Q=True, 
#                          degree = 0, ordering=["sens", "l"], l=2,
#                          cmap='grey',rot=0,
##                          savefig=True,
#                          )
#
## degree two
#SF = plot_coeffs_per_mode(ifiles,  label1 = "CC",
#                          ifiles2=ifiles2, label2="ellip CC", 
#                          mode='sc', model=["AD", "GLW", 'S20RTS+CRUST+BT'], 
#                          degree = 2, ordering="l",
#                          cmap='grey',rot=0, 
#                          l=2, #order=0, 
##                          savefig=True,
#                          )
#
## degree four
#SF = plot_coeffs_per_mode(ifiles,  label1 = "CC",
#                          ifiles2=ifiles2, label2="ellip CC", 
#                          mode='sc',  model=["AD", "GLW"],#,'S20RTS+CRUST+BT'], 
#                          degree = 4, ordering="l",
#                          cmap='grey',rot=0, 
#                          l=2,# order=0,
##                          savefig=True,
#                          )
#
## CC degrees
#SF = plot_coeffs_per_mode(ifiles,  label1 = "CC",
#                          mode='cc', model=['S20RTS','S20RTS+CRUST+BT'],
#                          degree=2, ordering="l",
#                          cmap='grey', 
#                          l=[0,2],
##                          savefig=True,
#                          )
#
##SF = plot_coeffs_per_mode(ifiles, 
##                          mode='cc', model=['S20RTS', "AD"],
##                          degree=4, ordering="l")
