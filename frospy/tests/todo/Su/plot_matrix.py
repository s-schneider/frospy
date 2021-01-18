#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 20:44:24 2018

@author: talavera
"""

#!/usr/bin/python

# #--------------
# from frospy.plot.nmplt import plot_Mmatrix  
# m = 'cstCatt'
# cwd = '/net/home/talavera/eejit/splitting/03s00/synthetics/IC/QCC'
# plot_Mmatrix(matrix_dir="%s/matrix.%s.dat"%(cwd,m), 
#              matrix_file='%s/omega.%s.dat'%(cwd,m), 
#              modes_dir='%s'%cwd, 
#              matrix_type='ascii',
# #             ylim=[230,500],
# #             vmin=-0.5, vmax=5,
# #             savefig=True,
#              )

# cwd = '/net/home/talavera/eejit/splitting/04s00/synthetics/IC/QCC'
# plot_Mmatrix(matrix_file="%s/matrix.dat"%(cwd), 
#              omega_file='%s/omega.dat'%(cwd), 
#              modes_dir='%s'%cwd, 
#              matrix_type='ascii',
# #             vmin=-0.5, vmax=5,
#              )


#-------------- My old version
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os
from frospy.util.read import read_omega_dat

#os.system('/data/talavera/src/bin_to_ascii_ww.py ww.dat')
#os.system('/data/talavera/src/bin_to_ascii_matrix.py') 
#os.system('/data/talavera/src/writeWW')
#
N = int(open('matrix.dat').readline().rstrip()) # matrix size from first row
M = np.loadtxt("matrix.dat", skiprows=1)
M = [ii for ii in zip(*M)]
#N = int(open('matrix.dat').readline().rstrip()) # matrix size from first row
#M = np.loadtxt("matrix.dat", skiprows=1); M = zip(*M)
M_re = np.array(M[0]).reshape((N,-1))#/1e-6
M_im = np.array(M[1]).reshape((N,-1))/1e-6
omega = read_omega_dat("omega.dat")
omega = [ii for ii in zip(*omega)]


modesin = open('modes.in').read().splitlines()
for i, mode in enumerate(modesin):
    modesin[i] = mode.split()

contrast = 1# to view low plaitudes in matrix


fig = plt.figure(figsize=(25,10))
ax = fig.add_subplot(1,3,1)
color = plt.cm.get_cmap("rainbow", int(modesin[0][0]))
l0 = 0; i = 0
for mode in modesin[1:int(modesin[0][0])+1]:
    n = int(mode[0])
    k = mode[1]
    l = int(mode[2])
    
    l1 = l0 + int(mode[2])*2 + 1
    ax.scatter(omega[0][l0:l1], omega[1][l0:l1], 
               label='${}_{%s}%s_{%s}$'%(n, k, l),
               marker='D', color=color(i), s = 50, alpha=0.95)
    l0 = l1; i = i + 1 
#ax.scatter(omega[0], omega[1], marker='^', color='g')
#ax.set_xlim(min(omega[0])*0.999, max(omega[0])*1.001)
#ax.set_xlim(9.83, 9.93)
#ax.set_ylim(785, 835)
aspect = (ax.get_xlim()[1] - ax.get_xlim()[0]) / \
        (ax.get_ylim()[1] - ax.get_ylim()[0])
ax.set_aspect(adjustable='box', aspect=aspect)
ax.set_xlabel("f(mHz)", fontsize=18)
ax.set_ylabel("Q", fontsize=18)
ax.legend(loc='upper center', bbox_to_anchor=(0.45, 1.15), 
          columnspacing=0,handletextpad=-0.5,
          ncol=len(modesin[1:]), fontsize=18, frameon=False)


ax = fig.add_subplot(1,3,2)
s = ax.imshow(M_re, cmap = plt.cm.get_cmap("jet", 25), aspect=1,
              vmin=min(np.min(M_re),np.min(M_im)), vmax=np.max(M_im)*contrast)
#              vmin=-15, vmax=4.5)
mf = ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2)) # for exponents greater than +-2
cb = fig.colorbar(s, orientation='horizontal', pad=0.01, format=mf)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.ax.set_xlabel('Re[M]($\mu$Hz)', fontsize=18)
ax.tick_params(axis='x',which='both', top='on', bottom='off',
               labelbottom='off', labeltop='on')
lsum= 0; m=[]; d = 0.5# d size of point in imshow
for mode in modesin[1:int(modesin[0][0])+1]:
    n = int(mode[0])
    k = mode[1]
    l = int(mode[2])
    lsum = lsum + (2*l + 1)
    m.extend(range(-l,l+1))
    ax.plot([lsum-d, lsum-d],[0-d,N-d], lw=1, c="w") # vertical
    ax.plot([0-d,N-d],[lsum-d, lsum-d],lw=1, c="w") # horizontal
ax.set_xticks(range(len(m)))
ax.set_xticklabels (m, fontsize=10)
ax.set_yticks(range(len(m)))
ax.set_yticklabels (m, fontsize=10)
ax.set_xlabel('$m$', fontsize=18)
ax.set_ylabel('$m^{\prime}$', rotation='vertical', fontsize=18)
ax.xaxis.set_label_position('top')
ax.set_xlim(-d, N-1+d)
ax.set_ylim(N-1+d, -d)


ax = fig.add_subplot(1,3,3)
s = ax.imshow(M_im, cmap = plt.cm.get_cmap("jet", 25), aspect=1,
              vmin=min(np.min(M_re),np.min(M_im)), vmax=np.max(M_im)*contrast)
#              vmin=-15, vmax=4.5)
mf = ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2)) # for exponents greater than +-2
cb = fig.colorbar(s, orientation='horizontal', pad=0.01, format=mf)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.ax.set_xlabel('Im[M]($\mu$Hz)', fontsize=18)
ax.tick_params(axis='x',which='both', top='on', bottom='off',
               labelbottom='off',labeltop='on')
lsum= 0; m=[]; d = 0.5# d size of point in imshow
for mode in modesin[1:int(modesin[0][0])+1]:
    n = int(mode[0])
    k = mode[1]
    l = int(mode[2])
    lsum = lsum + (2*l + 1)
    m.extend(range(-l,l+1))
    ax.plot([lsum-d, lsum-d],[0-d,N-d], lw=1, c="w") # vertical
    ax.plot([0-d,N-d],[lsum-d, lsum-d],lw=1, c="w") # horizontal
    ax.annotate('${}_{%s}%s_{%s}$'%(n, k, l), xy=(2, 1), 
                xytext=(N-1+d, lsum-l-d), 
                va='center', ha='left', fontsize=18)
ax.set_xticks(range(len(m)))
ax.set_xticklabels (m, fontsize=8)
ax.set_yticks(range(len(m)))
ax.set_yticklabels (m, fontsize=8)
ax.set_xlabel('$m$', fontsize=18)
ax.set_ylabel('$m^{\prime}$', rotation='vertical', fontsize=18)
ax.xaxis.set_label_position('top')
ax.set_xlim(-d, N-1+d)
ax.set_ylim(N-1+d, -d)
plt.show() 

#-----------------
#ata matrix
#M = np.loadtxt("ascii_atab.mat")
#print os.getcwd()
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#s = ax.imshow(M, cmap = plt.cm.get_cmap("jet", 25), aspect=1, vmin=-1e7, vmax=1e7)
#cb = fig.colorbar(s, orientation='horizontal', pad=0.01, format='%0.2f')
#tick_locator = ticker.MaxNLocator(nbins=4)
#cb.locator = tick_locator
#cb.update_ticks()
#ax.tick_params(axis='x',which='both', top='on', bottom='off',
#               labelbottom='off', labeltop='on')

