# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:07:45 2018

@author: sujania
"""

import numpy                               as np
import matplotlib.pyplot                   as plt
from mpl_toolkits.basemap                  import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
from frospy.plot.nmplt                       import beachball 

cmt_shallow = [
            '073095A',
            # '021796B', # not used in radial modes
            '032598B',
            '062301E', 
            # '092503C', # not used in radial modes
            '032805D',
            '111506F',
            '011307A',
            # '040107E', # not used in radial modes
            # '081507F', # not used in radial modes
            '022710A',
            '031111B',
            ]

cmt_deep = [
           # '120678A', # not used in radial modes
           '030994E',
           '060994A',
           '100494B',
           '061796A',
           '050306F',
           '070508A',
           '053015A', # new earthquake
           '052413A', # new earthquake
           ]

cmt_new = [
           '081918A', # not used in radial modes
           '090817A',
           '122516A',
           '100494B',
           '061796A',
           '050306F',
           '070508A',
           '053015A', # new earthquake
           '052413A', # new earthquake
           ]
cat_shallow = []
cat_deep = []
lon_0=0
for event in cmt_shallow: cat_shallow.append(read_std_cat(cmt_id=event))
for event in cmt_deep: cat_deep.append(read_std_cat(cmt_id=event))

fig = plt.figure(figsize=(8.5,5))
ax = fig.add_subplot(111)
m = Basemap(projection='kav7',lon_0=lon_0,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='lightblue')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='lightblue')
 
for cat in cat_shallow:
    cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
    fm  = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp]
    lat = cat[0].origins[0].latitude
    lon = cat[0].origins[0].longitude
    
    axins = inset_axes(ax, .1, .1)
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    x, y = m(lon, lat); h = .17/2.5; w = .075/2
    x = x/xmax - w/2.; y = y/ymax - h/2.
    ip = InsetPosition(ax, [x, y, w, h])
    axins.set_axes_locator(ip)
    beachball(fm, facecolor='r', linewidth=0.5, axis=axins)
ax.set_title('Shallow events (%s)'%len(cat_shallow),fontsize=25)
plt.tight_layout()
plt.show()
fig.savefig('cmts_shallow', dpi=350) 

fig = plt.figure(figsize=(8.5,5))
ax = fig.add_subplot(111)
m = Basemap(projection='kav7',lon_0=lon_0,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='lightblue')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='lightblue')
 
for cat in cat_deep:
    cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
    fm  = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp]
    lat = cat[0].origins[0].latitude
    lon = cat[0].origins[0].longitude
    
    axins = inset_axes(ax, .1, .1)
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    x, y = m(lon, lat); h = .17/2.5; w = .075/2
    x = x/xmax - w/2.; y = y/ymax - h/2.
    ip = InsetPosition(ax, [x, y, w, h])
    axins.set_axes_locator(ip)
    beachball(fm, facecolor='r', linewidth=0.5, axis=axins)
ax.set_title('Deep events (%s)'%len(cat_deep), fontsize=25)
plt.tight_layout()
plt.show()
fig.savefig('cmts_deep', dpi=350) 

cat_all = cat_shallow + cat_deep
fig = plt.figure(figsize=(8.5,5))
ax = fig.add_subplot(111)
m = Basemap(projection='kav7',lon_0=lon_0,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='lightblue')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='lightblue')
 
for cat in cat_all:
    cmt = cat[0].focal_mechanisms[0].moment_tensor.tensor
    fm  = [cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp]
    lat = cat[0].origins[0].latitude
    lon = cat[0].origins[0].longitude
    
    axins = inset_axes(ax, .1, .1)
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    x, y = m(lon, lat); h = .17/2.5; w = .075/2
    x = x/xmax - w/2.; y = y/ymax - h/2.
    ip = InsetPosition(ax, [x, y, w, h])
    axins.set_axes_locator(ip)
    beachball(fm, facecolor='r', linewidth=0.5, axis=axins)
ax.set_title('Events (%s)'%len(cat_all), fontsize=25)
plt.tight_layout()
plt.show()
fig.savefig('cmts_all', dpi=350) 