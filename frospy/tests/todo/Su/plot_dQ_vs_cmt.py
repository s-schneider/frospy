# -*- coding: utf-8 -*-
"""
Created on Sun May 20 13:49:37 2018

@author: sujania
"""

import numpy                               as np
import matplotlib.pyplot                   as plt
import matplotlib.cm                       as cm
from mpl_toolkits.basemap                  import Basemap
from mpl_toolkits.axes_grid1               import make_axes_locatable
from frospy.util.read                        import read_modes, read_std_cat
from frospy.preprocessing.spectrum           import get_mode
from frospy.splitting.coeff                  import calc_cst_indices   
from matplotlib.colors import LinearSegmentedColormap

cmt_id = [
       '060994A',
       '061796A',
       '030994E',
       '070508A',
       '022710A',
       '100494B',
       '031111B',
       '032805D',
       '111506F'
        ]

mode_name = "3S0"
modedata  = read_modes(); 
mode      = get_mode([mode_name], modedata) # reading mode data
Q0        = mode[mode_name]['Q']; 
f0        = mode[mode_name]['freq']*1e3
cat_3s0   = []

for event in cmt_id: cat_3s0.append(read_std_cat(cmt_id=event))

f = np.array([
                3272.5,
                3272.5,
                3272.6,
                3272.6,
                3269.9,
                3268.8,
                3271.1,
                3272.0,
                3271.4
                ])
Q = np.array([
                1155,
                1120,
                1336,
                1217,
                1214,
                4246,
                3106,
                1914,
                1897
                ])
#f0_inv = 3272.5 # Bolivia
#Q0_inv = 1136 # Bolivia
f0_inv = 3272.59 # HT
Q0_inv = 1217 # HT

df = f - f0_inv
dQ = Q - Q0_inv

fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(211)
m = Basemap(projection='kav7',lon_0=180,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(-90.,91.,60.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='white')

i = 0
for cat in cat_3s0:
    cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
    fm  = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp]
    lat = cat[0].origins[0].latitude
    lon = cat[0].origins[0].longitude
    x, y = m(lon, lat)
    scatter = m.scatter(x, y, c=df[i], vmin=-0.25, vmax=0.25,
                        s=500, marker="o", alpha=0.95, cmap=cm.coolwarm_r, zorder=10)
    i = i + 1

divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
fig.add_axes(cax)
clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
clb.ax.set_title('df', position=(0.5, -2), fontsize=18)
ax.set_title('${}_3S_0$', fontsize=30)

ax = fig.add_subplot(212)
m = Basemap(projection='kav7',lon_0=180,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(-90.,91.,60.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='white')

i = 0
for cat in cat_3s0:
    cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
    fm  = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp]
    lat = cat[0].origins[0].latitude
    lon = cat[0].origins[0].longitude
    x, y = m(lon, lat)
    scatter = m.scatter(x, y, c=dQ[i], vmin=-400, vmax=400,
                        s=500, marker="o", alpha=0.95, cmap=cm.coolwarm_r, zorder=10)
    i = i + 1

divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
fig.add_axes(cax)
clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
clb.ax.set_title('dQ', position=(0.5, -2), fontsize=18)
plt.tight_layout()
plt.show()    
#fig.savefig('3s0_df_dQ_wrt_HT_deepevents', dpi=350) 