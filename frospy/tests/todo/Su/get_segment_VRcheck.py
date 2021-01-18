#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:55:37 2018

@author: talavera
"""

import sys
import os, fnmatch
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print "USAGE python bin_to_ascii.py damp"
    sys.exit()

cwd = os.getcwd() # pwd
#event = sys.argv[1] # Event name
damp = float(sys.argv[1])
it_max = len(fnmatch.filter(os.listdir('%s/d%.0f' % (cwd, damp)), 'Iter_*'))

event_file = '%s/d%.0f/config/events_list.dat' %(cwd, damp)
event_name = os.path.basename(cwd)
event_list = open(event_file).read().splitlines()[1:]

events = [[],[]]
segments = [[],[]]

for e in event_list:
    event_file = '%s/d%.0f/Iter_%d/synseis/%s_tot_misfit.dat' % (cwd, damp, it_max, e)
    segment_file = '%s/d%.0f/Iter_%d/synseis/%s_seg_misfit.dat' % (cwd, damp, it_max, e)
    events[0].append(float(open(event_file).readlines()[0].split()[1]))
    events[1].append(float(open(event_file).readlines()[1].split()[1]))
    segments[0].append(list(map(str.strip, open(segment_file).readlines()[0::6]))) 
    segments[1].append(map(float,str(open(segment_file).readlines()[4::6]).split()[1:][1::4]))

print event_name, '\t', events[1]

fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(1,1,1)
ax.barh(range(2*len(event_list))[1::2],events[1], 
        color='b', height=0.58, align='center', label='stations')
#ax.set_xlabel('# stations', fontsize=10)
ax.set_xlim(-max(events[1])*1.1, max(events[1])*1.1)
ax.set_xticks(np.arange(0, max(events[1])*1.1, 25))
ax.set_xlabel('# stations', horizontalalignment='right', position=(0.8,25), fontsize=10)

ax = ax.twiny()
ax.barh(range(2*len(event_list))[0::2],events[0],
        color='r', height=0.69, align='edge', label='Var. Red')
ax.set_yticks(np.linspace(0.5, 2*len(event_list)+0.5,2*len(event_list)+1)[0::2])
ax.set_yticklabels (event_list, rotation='horizontal', fontsize=10)
ax.set_xlim(-max(events[0])*1.1, max(events[0])*1.1)
ax.set_xlabel('Variance Reduction', fontsize=10)
vals = ax.get_xticks(); ax.set_xticklabels(['{:3.0f}%'.format(x) for x in vals])
fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.suptitle('%s with damp=%s'%(event_name, 1./damp))
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.show()

i = 0; key = ['nps']
while key != 'quit':
#for i in range(6):
    if len(segments[0][i]) > 100: y = 25
    elif len(segments[0][i]) > 50 and len(segments[0][i]) < 100: y = 15
    elif len(segments[0][i]) < 50 and len(segments[0][i]) > 10: y = 10
    elif len(segments[0][i]) < 10: y = 5
    
    #print len(segments[0][i]) , y
    
    fig = plt.figure(figsize=(7,y))
    ax = fig.add_subplot(1,1,1)
    ax.barh(range(len(segments[0][i])),segments[1][i],
            color='b', height=0.69, align='center', label='Var. Red')
    ax.set_yticks(range(len(segments[0][i])))
    ax.set_yticklabels (segments[0][i], rotation='horizontal', fontsize=10)
    ax.set_xlim(0, 110)
    vals = ax.get_xticks()
    ax.set_xticklabels(['{:3.0f}%'.format(x) for x in vals])
    ax.plot(np.ones(len(segments[0][i])+2)*events[0][i], range(-1,len(segments[0][i])+1),
            color='r', linestyle='--', linewidth=3)
    ax.set_ylim(-1, len(segments[0][i]))
    fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
    ax.set_title(r'%s' % event_list[i], fontsize=15)
    plt.show()    
    
    print 'Event %s: %d/%d with %d stations' %(event_list[i], i+1, len(event_list), len(segments[0][i]))
    key = raw_input('Please an action (nsp, bts, psp, quit): ')
    if i >= len(event_list)-1 : i = 0
    elif key == 'bts': i = 0
    elif key == 'next' or key == 'nsp' or key == "": i += 1
    elif key == 'quit' or key == 'exit' or key == 'q': break
    elif key == 'psp':
        if i != 0:
            i -= 1
        else:
            i = 0
            print('Already at the beginning.\n')