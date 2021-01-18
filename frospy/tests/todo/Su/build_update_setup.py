#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:33:33 2019

@author: talavera
"""
import os
import pickle
from os.path import join, basename
from frospy.core.setup.settings import update_setup_pickle
from frospy.util.read import read_pickle

# location of setups: find . -name "setup.pickle” 
net = "/net/home/talavera/eejit"
quanta = "/quanta1/home/talavera"
iloc  = "splitting/02s00-07s02"
ifile = "allevents/cst+d00/cst_02s00"
setup_files = [join(net, iloc, ifile, "setup.pickle")]
for s in setup_files:
#    setup = pickle.load(open(s, "rb" ))
    setup = read_pickle(s)
    model = setup.model
    if model.startswith("//nfs"):
        print("changing mzero")
        mzero = join(net, iloc, ifile, "mzero.dat")
        inversion_out = join(net, iloc, ifile, "inversion_out")
        os.system('cp %s %s' % (mzero, inversion_out))
        
        setup.model  = "PREM"
        setup.rundir = join(quanta, iloc, ifile)
        setup.bindir = join(quanta, "bin")
        setup.nmpydir= join(quanta, "codes/nmPy")
        for key, value in setup.segmentsdir.iteritems():
            segname = basename(value)
            setup.segmentsdir[key] = join(quanta, iloc, "segments", segname)
        path = os.path.dirname(s)
        setup = update_setup_pickle(setup)
        setup.write(path)#, overwrite=True)
    else:
        setup.model  = "PREM"
        path = os.path.dirname(s)
        setup = update_setup_pickle(setup)
        setup.write(path)#, overwrite=True)
    

#from frospy.core.segment import read
#
#
#sdir = '/net/home/talavera/eejit/etc/tests/segments/111401B.dat'
#
#seg = read(sdir)
#
#seg.sort()
#seg.picks