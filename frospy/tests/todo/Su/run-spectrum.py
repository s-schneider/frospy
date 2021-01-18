#!/usr/bin/python

# Allows to run spectrum and saves data in ASCII
from frospy.preprocessing.spectrum import spectrum
from frospy.preprocessing.spectrum import taper
from frospy.preprocessing.spectrum import printw
import sys

event = sys.argv[1] # Event name

spectrum(data='%s.ahx'%event, syn='%s.ahx.syn.fil'%event, show_modes=True)

