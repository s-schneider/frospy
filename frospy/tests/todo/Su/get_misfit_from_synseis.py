#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 19:19:58 2018

@author: talavera
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from frospy.util.read import read_std_cat
from frospy.util.read import read_modes
from frospy.preprocessing.spectrum import get_mode

#inv   = read_std_inv()
event ='060994A'
#cat   = read_std_cat(event)
cwd   = os.getcwd()

mode_name = os.path.basename(cwd)
modedata  = read_modes(); 
mode      = get_mode([mode_name], modedata) # reading mode data

for name,dict_ in mode.items(): mode_name = name
Q0        = mode[mode_name]['Q']; 
f0        = mode[mode_name]['freq']*1e3

# data
re = []; i=0; lb=0; stations=[];tw=[] 
with open('inv-self/%s_obs_r.dat'%event,'rb') as obs_r:
    for line in obs_r:
        lines = line.split('\t')
        if i == 0: stations.append(lines[0])
        if i == 2: tw = lines[0].split()
        if i == 3: lb = int(lines[0])
        if i > 5:
            re.append(complex(lines[0]))
        i = i + 1
        if i == lb+6: i = 0
        
re = [re[lb*j:lb*(j+1)] for j in range(len(re)/lb)]

im = []; i=0; lb=0
with open('inv-self/%s_obs_i.dat'%event,'rb') as obs_i:
    for line in obs_i:
        lines = line.split('\t')
        if i == 3: lb = int(lines[0])
        if i > 5:
            im.append(float(lines[0])* 1j)
        i = i + 1
        if i == lb+6: i = 0
im = [im[lb*j:lb*(j+1)] for j in range(len(im)/lb)]

re_im_obs = []
for i in range(len(im)):
    re_im_obs.append(np.array(im[i]) + np.array(re[i]))

# syn
#inversion =  ['HT', 'QM1', 'DE', 'SS', 'inv-self', 'Qtest']
#inversion =  ['HT', 'QM1', 'DE', 'SS', 'inv-self', 'inv-cross']
inversion =  ['HT', 'QM1', 'SS', 'inv-self', 'inv-cross']
syn=[]; misfit=[]; amp_ratio = []; tot_VR=[[],[]]; tot_misfit=[[],[]]
for model in inversion:
    re = []; i=0; lb=0
    with open('%s/%s_syn_r.dat'%(model, event),'rb') as syn_r:
        for line in syn_r:
            lines = line.split('\t')
            if i == 3: lb = int(lines[0])
            if i > 5:
                re.append(complex(lines[0]))
            i = i + 1
            if i == lb+6: i = 0
            
    re = [re[lb*j:lb*(j+1)] for j in range(len(re)/lb)]
    
    im = []; i=0; lb=0
    with open('%s/%s_syn_i.dat'%(model, event),'rb') as syn_i:
        for line in syn_i:
            lines = line.split('\t')
            if i == 3: lb = int(lines[0])
            if i > 5:
                im.append(float(lines[0])* 1j)
            i = i + 1
            if i == lb+6: i = 0
            
    im = [im[lb*j:lb*(j+1)] for j in range(len(im)/lb)]

    re_im_syn = []
    for i in range(len(im)):
        re_im_syn.append(np.array(im[i]) + np.array(re[i]))

    
    mf = [[],[]]; ratio=[]; r_i_obs=[]; r_i_syn=[]
    for i in range(len(re_im_obs)):
        #misfit
        d2=sum(abs(re_im_obs[i])**2.)
        rmf = sum((abs(re_im_obs[i]) - abs(re_im_syn[i]))**2.) / d2
        cmf = sum(abs(re_im_obs[i]- re_im_syn[i])**2.) / d2
        #mf[0].append((1-rmf)*100); mf[1].append((1-cmf)*100) # VR
        mf[0].append(rmf); mf[1].append(cmf) # misfit
        r_i_obs.extend(list(re_im_obs[i]))
        r_i_syn.extend(list(re_im_syn[i]))
    
        # amplitude ratio=data/syn
        ratio.append(max(abs(re_im_obs[i]))/max(abs(re_im_syn[i])))

    syn.append(re_im_syn)
    misfit.append(mf)
    amp_ratio.append(ratio)
    #tot_VR[0].append(float(open('%s/Total_VR.out'%model,'rb').readlines()[2].split(':')[1]))
    #tot_VR[1].append(float(open('%s/Total_VR.out'%model,'rb').readlines()[1].split(':')[1]))
    
    r_i_obs=np.array(r_i_obs); r_i_syn=np.array(r_i_syn)
    d2=sum(abs(r_i_obs)**2.)
    rmf = sum((abs(r_i_obs) - abs(r_i_syn))**2.) / d2
    cmf = sum(abs(r_i_obs- r_i_syn)**2.) / d2
    tot_misfit[0].append(rmf); tot_misfit[1].append(cmf)

 
#plot complex misfit
#color = ['b', 'g', 'c', 'plum', 'm', 'r']; 
color = ['b', 'g', 'c', 'plum', 'r']; 
#color = ['b', 'g', 'c', 'r']; 
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(3,1,1)
i=0
for m in misfit:
    ax.plot(m[1], color[i], label=r'%s=%s'%(inversion[i], round(tot_misfit[1][i],3)))#round(tot_VR[1][i],1)))
    i = i + 1
ax.set_title('$%s$ with %s stations from %s, tw=%s-%s'
             %(mode[mode_name]['name'], len(stations), event, float(tw[0]), float(tw[1])), fontsize=18)
#ax.set_ylabel('Complex VR',fontsize=12)
ax.set_ylabel('Complex misfit',fontsize=12)
ax.set_xticks(range(len(stations)))
ax.set_xticklabels([])
ylim=ax.get_ylim()
plt.legend(loc=4)

#plot amplitude misfit
ax = fig.add_subplot(3,1,2)
i=0
for m in misfit:
    ax.plot(m[0], color[i], label=r'%s=%s'%(inversion[i], round(tot_misfit[0][i],3)))#round(tot_VR[0][i],1)))
    i = i + 1
#ax.set_ylabel('Amplitude VR', fontsize=12)
ax.set_ylabel('Amplitude misfit', fontsize=12)
ax.set_xticks(range(len(stations)))
ax.set_xticklabels([])
#ax.set_ylim(ylim)
plt.legend(loc=4)

#plot amplitude ratio
ax = fig.add_subplot(3,1,3)
i=0
for r in amp_ratio:
    ax.plot(r, color[i], label=r'%s=%s'%(inversion[i], round(np.mean(r),2)))
    i = i + 1
ax.plot(np.ones(len(stations)), '--k', linewidth=0.75)
ax.set_ylabel('data/syn ratio', fontsize=12)
ax.set_xticks(range(len(stations)))
ax.set_xticklabels (stations, rotation='vertical', fontsize=10)
ax.invert_yaxis()
plt.tight_layout()
plt.legend(loc=4)
plt.show()

