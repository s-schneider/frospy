#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:49:36 2018

@author: talavera
"""

# Allows to run spectrum and saves data in ASCII
from frospy.preprocessing.spectrum import spectrum
import sys
import os, fnmatch

print 'the script should be run in the folder contaning the damped iterations'
if len(sys.argv) < 3:
    print "USAGE python bin_to_ascii.py event mode(i.e 01s00)"
    sys.exit()

#event='060994A'
#mode='01s00'
event = sys.argv[1] # Event name
mode = sys.argv[2] # mode name
cwd = os.getcwd() # pwd
#cwd = '%s/%s' %(os.getcwd(), event) # pwd

data = '/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx'%event
s20rts = '/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx.syn'%event
segment = '/net/home/talavera/radial-inversion/%s/%s/%s/%s.segments' %(mode, mode, event,event)
damp = fnmatch.filter(os.listdir(cwd), 'd1*') + fnmatch.filter(os.listdir(cwd), 'd0*')
it_max = len(fnmatch.filter(os.listdir(cwd+'/'+damp[0]), 'Iter_*'))
fw = map(float, open('%s'%segment).readlines()[1].replace('\t',' ').replace('\n',' ').split())
tw = map(float, open('%s'%segment).readlines()[2].replace('\t',' ').replace('\n',' ').split())

syn = [s20rts]; label = ['0-7.5mHz']

#previous models
for m in ['HT', 'QM1', 'DE']:
   syn.append('%s/%s/Iter_1/synseis/%s.ahx.syn' % (cwd, m, event))
   label.append('%s, model=%s' % (mode, m))
   
# damped synthetics
for d in sorted(damp):
   syn.append('%s/%s/Iter_%d/synseis/%s.ahx.syn' % (cwd, d, it_max, event))
   label.append('%s, damp=%s' % (mode, str(1/float(d[1:]))))
   
# Only self
spectrum(data=data, syn=syn, syn_label=label, segment_file=segment, 
         fw=[fw[0]*0.98, fw[1]*1.02], tw=tw, show_modes=True)

#os.remove('%s/%s.segments' %(cwd, event))