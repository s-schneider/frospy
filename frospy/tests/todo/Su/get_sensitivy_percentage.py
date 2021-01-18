#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:45:59 2019

@author: talavera
"""

import numpy as np
import glob
from frospy.core.modes import read as read_modes
from frospy.core.splittingfunc.plot import sens_kernel
import matplotlib.pyplot as plt

chome = '/net/home/talavera/forward/splitting/dat-kernel'
modelist = [
#'002s003',
#'003s001',
#'003s002',
#'005s002',
#'006s003',
#'008s001',
#'008s002',
#'008s005',
#'009s002',
#'009s003',
#'009s004',
#'010s002',
#'011s002',
#'011s004',
#'011s005',
#'011s006',
#'013s001',
#'013s002',
#'013s003',
#'013s006',
'014s004',
#'015s003',
#'015s004',
#'016s005',
#'016s006',
#'016s007',
#'017s001',
#'017s008',
#'018s003',
#'018s004',
#'018s006',
#'020s001',
#'020s005',
#'021s006',
#'021s007',
#'022s001',
#'023s004',
#'023s005',
#'025s001',
#'025s002',
#'027s002',
]

for mode in modelist:
    # vp, rho, vs
#    mode = '027s002'
    ifiles = glob.glob('%s/*kernel*%s*'%(chome,mode))
    cfile = []
    
    r_TZ = 5700
    r_CMB = 3400
    r_ICB = 1220
    vs = 0
    vp = 0
    
    allmodes = read_modes()
    m = read_modes(modenames=mode)[0]
    print(m.name)
    
    for f in ifiles:
        cfile.append(open(f).readlines())
        c = open(f).read().splitlines()
        c = [line.strip().split() for line in c]
        c = zip(*c)
        r = [float(e) for e in c[0]]
        sens = [(float(e)) for e in c[1]]
        kind = f.split('/')[-1].split('-')[0]
        
        rr_CMB = next(x for x in r if x > r_CMB)
        i_CMB = r.index(rr_CMB)
        
        rr_ICB = next(x for x in r if x > r_ICB)
        i_ICB = r.index(rr_ICB)
        
    #    total  = sum(sens)
    #    mantle = sum(sens[i_CMB::])/total
    #    outerc = sum(sens[i_ICB:i_CMB])/total
    #    innerc = sum(sens[0:i_ICB])/total
        total = np.trapz(sens,x=r,dx=0.1)
        mantle = np.trapz(sens[i_CMB::],x=r[i_CMB::],dx=0.1)/total
        outerc = np.trapz(sens[i_ICB:i_CMB],x=r[i_ICB:i_CMB],dx=0.1)/total
        innerc = np.trapz(sens[0:i_ICB],x=r[0:i_ICB],dx=0.1)/total
#        if abs(mantle) > 1:
#            mantle = 0
#        if abs(outerc) > 1:
#            outerc = 0
#        if abs(innerc) > 1:
#            innerc = 0
        
        if kind == 'vs':
            vs = innerc*total
        elif kind == 'vp':
            vp = innerc*total
        
        if kind == 'vp' or kind == 'vs':
            print("%s:\t mantle: %s\t oc: %s\t ic: %s\t total:%s"
                   %(kind, 
                     round(mantle,3)*100, 
                     round(outerc,3)*100, 
                     round(innerc,3)*100,
                     round(mantle+outerc+innerc,3)*100))

    #    print(np.trapz(sens[i_CMB::],x=r[i_CMB::]), 
    #          np.trapz(sens[i_ICB:i_CMB],x=r[i_ICB:i_CMB]), 
    #          np.trapz(sens[0:i_ICB],x=r[0:i_ICB]),
    #          np.trapz(sens,x=r,))
    #    fig = plt.figure()
    #    ax = fig.add_subplot(1,1,1)
    #    ax.plot(sens[i_CMB::],r[i_CMB::])
    #    ax.plot(sens[i_ICB:i_CMB],r[i_ICB:i_CMB])
    #    ax.plot(sens[0:i_ICB],r[0:i_ICB])
    #    ax.set_title(kind)
    if vs > vp:
        print('%s: Vs mode'%m.name)    
    print('\n')
    sens_kernel(m, title=mode)