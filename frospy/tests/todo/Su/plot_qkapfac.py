#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:04:50 2019

@author: talavera
"""
import os
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from os.path import join, basename, splitext

greens = plt.cm.get_cmap("Greens", 15)
blues = plt.cm.get_cmap("Blues", 15)

color = [blues(10), blues(7), blues(2), #blues(1), 
#         greens(2), greens(4), greens(6), greens(8), greens(10),greens(12),]
         greens(8), greens(10),greens(12),]
#         greens(2), greens(4), greens(6), greens(8), ]

dhome = "/net/home/talavera/eejit/splitting/"
mdir = ["%s/synthetics/Qkappa/fQ_qkapfac_IC.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_UM.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_lith.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_LM.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_OC.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_RT.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_QL6.dat",
        "%s/synthetics/Qkappa/fQ_qkapfac_QL6test.dat",
        ]

mode = ['00s00','01s00','02s00-07s02',
        '03s00','04s00','05s00','06s00',
        '07s00','08s00','09s00','11s00-27s02',]

ifiles_IC = [glob.glob(join(dhome,mdir[0]%m))[0] for m in mode]
ifiles_UM = [glob.glob(join(dhome,mdir[1]%m))[0] for m in mode]
ifiles_lith = [glob.glob(join(dhome,mdir[2]%m))[0] for m in mode]
ifiles_LM = [glob.glob(join(dhome,mdir[3]%m))[0] for m in mode]
ifiles_OC = [glob.glob(join(dhome,mdir[4]%m))[0] for m in mode]
ifiles_RT = [glob.glob(join(dhome,mdir[5]%m))[0] for m in mode]
ifiles_QL6 = [glob.glob(join(dhome,mdir[6]%m))[0] for m in mode]
ifiles_QL6test =[glob.glob(join(dhome,mdir[7]%m))[0] for m in mode]

Q_RM =  [5327,1499,1242,1083, 969, 921, 913,881, 852,840,832]
#Q_CC =  [5983,1856,1788,1127,1181,1006,1083,905,1028,937,985]
#err_u = [40, 144, 113, 17,  85,  19,  103, 121,75,  84, 58]
#err_l = [127, 143, 122, 13,  34,  8,   43,  7,  33,  45, 49]
Q_CC =  [5983,1856,1788,1228,1181,1006,1083,905,1028,937,985]
err_u = [40, 144, 113, 318,  85,  19,  103, 121,75,  84, 58]
err_l = [127, 143, 122, 35,  34,  8,   43,  7,  33,  45, 49]


#3s0 SC
#f=3272.40 -0.04/ +0.06
#Q=Q=1222.47 -22.07/ +61.88

#3s0 GC
#f=3271.93 -0.03/ +0.16
#Q=1228.17 -34.66/ +73.82


# Resovsky & Trampert,2005 /QL6
fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0
for f,Q,Q0,u,l in zip(ifiles_RT,Q_CC,Q_RM,err_u,err_l): #RT
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    if i == 0:
        ax.scatter(i,fQ[1]-Q0, marker='s',color='gray',
                   label="RT05: $Q_\kappa$=1e10@LM; $Q_\kappa^{max}$ elsewhere")
        ax.errorbar(i,Q-Q0, color="r", yerr=eval, marker='o', capsize=3,
                    label="CC", )
    ax.scatter(i,fQ[1]-Q0, marker='s',color='gray')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1

i = 0
for f,Q0 in zip(ifiles_QL6,Q_RM): #Ql6
    fQ = np.loadtxt(f)
    if i == 0:
        ax.scatter(i,fQ[1]-Q0, marker='s',color='k',
                   label="QL6: $Q_\kappa$=943@UM; $Q_\kappa$=1e10 elsewhere")
    ax.scatter(i,fQ[1]-Q0, marker='s',color='k')
    i = i + 1
    
i = 0
for f,Q0 in zip(ifiles_QL6test,Q_RM): #Ql6 test
    fQ = np.loadtxt(f)
    if i == 0:
        ax.scatter(i,fQ[1]-Q0, marker='s',color='b',
                   label="QL6: $Q_\kappa$=1328@UM")
    ax.scatter(i,fQ[1]-Q0, marker='s',color='b')
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Other models")
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
#ax.set_ylim(-1000,1200)
plt.legend(ncol=1,loc="upper right")



# Inner core
Qkap_IC = 1328.
#Qkap = [1000,1328.,2000,4427,10000,15000,30000,60000]
Qkap = [300,1000,1328.,3500,27000,60000]
#Qkap.extend(range(2000,60000,10000))
fac= [1./(Q/Qkap_IC) for Q in Qkap]
#for Q,f in zip(Qkap,fac):
##    print("IC with Qkappa=%s"%Qkap_IC,Q,f)
#    print(f)
##print(len(fac))

fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0
for f,Q,Q0,u,l in zip(ifiles_IC,Q_CC,Q_RM,err_u,err_l):
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    for p,c,Qk in zip(fQ,color,Qkap):
#        print(p)
        if i == 0:
            ax.scatter(i,p[1]-Q0, marker='s',
                        color=c, label="IC: %.2f, %s"%(1/p[2],Qk))
        ax.scatter(i,p[1]-Q0, marker='s',color=c)
    if i == 0:
        ax.errorbar(i,Q-Q0, color="r", label="CC", yerr=eval, marker='o')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Varying Inner Core $Q_\kappa$=%d"%Qkap_IC)
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
ax.set_ylim(-1000,1200)
plt.legend(ncol=2,loc="lower right")



# Upper mantle
Qkap_UM = 57822.
Qkap = [300,1000,1328.,3500,27000,60000]
#Qkap.extend(range(2000,60000,10000))
fac= [1./(Q/Qkap_UM) for Q in Qkap]
#for Q,f in zip(Qkap,fac):
##    print("UM with Qkappa=%s"%Qkap_UM,Q,f) 
#    print(f) 
##print(len(fac))

fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0
for f,Q,Q0,u,l in zip(ifiles_UM,Q_CC,Q_RM,err_u,err_l):
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    for p,c,Qk in zip(fQ,color,Qkap):
        if i == 0:
            ax.scatter(i,p[1]-Q0, marker='s',
                        color=c, label="UM: %.3f, %s"%(1/p[2],Qk))
        ax.scatter(i,p[1]-Q0, marker='s',color=c)
    if i == 0:
        ax.errorbar(i,Q-Q0, color="r", label="CC", yerr=eval, marker='o')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Varying Upper Mantle $Q_\kappa$=%d with IC $Q_\kappa$=%d"
             %(Qkap_UM,Qkap_UM))
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
ax.set_ylim(-1000,1200)
plt.legend(ncol=2,loc="lower right")



# Lower Mantle
#color = [blues(10), blues(7), blues(5), #blues(1), 
#         greens(4), greens(6), greens(8), greens(10),greens(12),]
##         greens(2), greens(4), greens(6), greens(8), ]

Qkap_LM = 57822.
Qkap = [300,1000,1328.,3500,27000,60000]
#Qkap = [300,1000,2000,4000,10000,15000,30000,60000]
#Qkap.extend(range(2000,60000,10000))
fac= [1./(Q/Qkap_LM) for Q in Qkap]
#for Q,f in zip(Qkap,fac):
##    print("Lower Mantle with Qkappa=%s"%Qkap_LM,Q,f) 
#    print(f) 
##print(len(fac))

fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0

for f,Q,Q0,u,l in zip(ifiles_LM,Q_CC,Q_RM,err_u,err_l):
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    for p,c,Qk in zip(fQ,color,Qkap):
#        print(p)
        if i == 0:
            ax.scatter(i,p[1]-Q0, marker='s',
                        color=c, label="LM: %.3f, %s"%(1/p[2],Qk))
        ax.scatter(i,p[1]-Q0, marker='s',color=c)
    if i == 0:
        ax.errorbar(i,Q-Q0, color="r", label="CC", yerr=eval, marker='o')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Varying Lower Mantle $Q_\kappa$=%d with IC $Q_\kappa$=%d"
             %(Qkap_LM,Qkap_LM))
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
ax.set_ylim(-1000,1200)
plt.legend(ncol=2,loc="lower right")



# Outer Core
#color = [blues(10), blues(7), blues(5), #blues(1), 
#         greens(4), greens(6), greens(8), greens(10),greens(12),]
##         greens(2), greens(4), greens(6), greens(8), ]

Qkap_OC = 57822.
Qkap = [300,1000,1328.,3500,27000,60000]
#Qkap = [300,1000,2000,4000,10000,15000,30000,60000]
#Qkap = [30,300,1000,2000,4000,10000,30000,60000]
#Qkap.extend(range(2000,60000,10000))
fac= [1./(Q/Qkap_OC) for Q in Qkap]
#for Q,f in zip(Qkap,fac):
##    print("Outer Core with Qkappa=%s"%Qkap_OC,Q,f) 
#    print(f) 
##print(len(fac))

fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0
for f,Q,Q0,u,l in zip(ifiles_OC,Q_CC,Q_RM,err_u,err_l):
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    for p,c,Qk in zip(fQ,color,Qkap):
#        print(p)
#        print(i,p[1],Q0,p[1]-Q0)
        if i == 0:
            ax.scatter(i,p[1]-Q0, marker='s',
                        color=c, label="OC: %.2f, %s"%(1/p[2],Qk))
        ax.scatter(i,p[1]-Q0, marker='s',color=c)
    if i == 0:
        ax.errorbar(i,Q-Q0, color="r", label="CC", yerr=eval, marker='o')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Varying Outer Core $Q_\kappa$=%d with IC $Q_\kappa$=%d"
             %(Qkap_OC,Qkap_OC))
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
ax.set_ylim(-1000,1200)
plt.legend(ncol=2,loc="lower right")




# Lithosphere
Qkap_lith = 57822.
Qkap = [300,1000,1328.,3500,27000,60000]
#Qkap = [300,1000,2000,4000,10000,15000,30000,60000]
#Qkap = [30,300,1000,2000,4000,10000,30000,60000]
#Qkap.extend(range(2000,60000,10000))
fac= [1./(Q/Qkap_lith) for Q in Qkap]
#for Q,f in zip(Qkap,fac):
##    print("Lithosphere with Qkappa=%s"%Qkap_lith,Q,f) 
#    print(f) 
##print(len(fac))

fig = plt.figure()
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1,1,1)
i = 0
for f,Q,Q0,u,l in zip(ifiles_lith,Q_CC,Q_RM,err_u,err_l):
    fQ = np.loadtxt(f)
    eval = np.array([l, u]).reshape((2, 1))
    for p,c,Qk in zip(fQ,color,Qkap):
#        print(p)
        if i == 0:
            ax.scatter(i,p[1]-Q0, marker='s',
                        color=c, label="Lith: %.2f, %s"%(1/p[2],Qk))
        ax.scatter(i,p[1]-Q0, marker='s',color=c)
    if i == 0:
        ax.errorbar(i,Q-Q0, color="r", label="CC", yerr=eval, marker='o')
    ax.errorbar(i,Q-Q0, color="r", yerr=eval,marker='o',capsize=3)
    i = i + 1
ax.axhline(y=0, color="gray", zorder=0)
ax.set_title("Varying Lythosphere $Q_\kappa$=%d with IC $Q_\kappa$=%d"
             %(Qkap_lith,Qkap_lith))
modename = ['${}_{%s}S_0$'%int(ll) for ll in range(0,10)]
modename.append('${}_{11}S_0$')
ax.set_xticks(np.linspace(0, len(modename), len(modename)+1))
ax.set_xticklabels(modename, rotation='horizontal', fontsize=15)
ax.set_xlim(-0.5,10.5)
ax.set_ylim(-1000,1200)
plt.legend(ncol=2,loc="upper right")