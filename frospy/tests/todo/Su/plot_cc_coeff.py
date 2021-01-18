# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 11:15:01 2018

@author: sujania
"""

import numpy as np
import matplotlib.pyplot as plt

# 05s00-13s02
cc_inv = np.array([-1.024446e+01, -1.111985e-01, -4.172265e+00, 7.545845e+00, -2.760083e+00])
cc_wh  = np.array([-23.249747946970412, 0., 0., 0., 0.])
cc_trp = np.array([-7.8685474666240403, 0., 0., 0., 0.])
cc_bgt = np.array([-23.243773926407343, 0., 0., 0., 0.])
cc_rwz = np.array([-19.445862735380029, 0., 0., 0., 0.])
cc_s20rts = np.array([-1.5207712954483830, -0.79245287657720009, 
             -1.8023340224734921, 2.5102791743920361, -3.7194176036023707])
cc_trmp_s20rts = cc_trp + cc_s20rts
cc_bgt_s20rts  = cc_bgt + cc_s20rts

fig = plt.figure(figsize=(10,3.5))
ax = fig.add_subplot(1,1,1)
ax.plot(cc_inv, '-^', color='r', label='${}_5S_0-{}_{13}S_2$', linewidth=2.5)
ax.plot(cc_bgt_s20rts, '-^', color='b', label='IC Beghein-Trampert + S20RTS', linewidth=2)
#ax.plot(cc_s20rts, '-^', color='b', label='S20RTS')
#ax.plot(cc_wh, '-^', color='g', label='IC Woodhouse')
#ax.plot(cc_bgt, '-^', color='c', label='IC Beghein-Trampert')
#ax.plot(cc_trp, '-^', color='plum', label='IC Tromp')
#ax.plot(cc_rwz, '-^', color='m', label='IC Romanowicz')
ax.set_xticks(range(len(cc_inv)))
plt.axhline(y=0, color='darkgray', linestyle='-',  linewidth=5, zorder=0)
plt.text(2, 1.5, 'PREM', fontsize=15, backgroundcolor='w', va='bottom', ha='center', zorder=0)
ax.tick_params(axis = 'both', which = 'major', labelsize = 15, zorder=0)
cst = ['$Re[c_{20}]$', '$Re[c_{21}]$', '$Im[c_{21}]$', '$Re[c_{22}]$', '$Im[c_{22}]$']
ax.set_xticks(range(len(cst)))
ax.set_xticklabels(cst, rotation='horizontal', fontsize=15)

# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
# Put a legend below current axis
ax.legend(loc='lower right',
          fancybox=True, shadow=True, ncol=2, fontsize=16)
#plt.title('Cross-coupling coefficients')
fig.savefig('cc_coeffs_05s00_13s02', dpi=350) 
