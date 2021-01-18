#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 14:45:12 2018

@author: talavera
"""
#%matplotlib inline 
#%matplotlib auto #plot outside of line
import os, sys
import numpy                               as np
import matplotlib.pyplot                   as plt
import matplotlib.cm                       as cm
from mpl_toolkits.basemap                  import Basemap
from mpl_toolkits.axes_grid1               import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
from frospy.plot.nmplt                       import beachball
from frospy.util.read                        import read_modes, read_std_cat
from frospy.preprocessing.spectrum           import get_mode
from frospy.splitting.coeff                  import calc_cst_indices
#from frospy.util.base                        import sc_degrees, cc_degrees

#if len(sys.argv) < 3:
#    print "USAGE python bin_to_ascii.py damp"
#    sys.exit()

cwd       = os.getcwd() # pwd
#damp      = float(sys.argv[1])
#mode_name = float(sys.argv[2])
damp      = 1000
mode_name = "1S0"
modedata  = read_modes(); 
mode      = get_mode([mode_name], modedata) # reading mode data
Q0        = mode[mode_name]['Q']; 
f0        = mode[mode_name]['freq']*1e3
cat       = []

cmt_id = open('%s/allevents/d%s/config/events_list.dat' %(cwd, damp)).read().replace('\n',' ').split()[1:]
for event in cmt_id: cat.append(read_std_cat(cmt_id=event))

# getting cst indices automatically
Q_file    = "%s/allevents/d%s/Iter_10/cst_solver/new_model.damp_%.2e.cst" %(cwd, damp, damp)
cst_ind   = calc_cst_indices('%s/allevents/d%s/config' %(cwd, damp))
cst_ind   = zip(*cst_ind); cst_ind[0] = np.array(cst_ind[0]) - 1 # setting indexes as python
mode_ii   = [i for i, x in enumerate(cst_ind[1]) if x == mode_name]
mode_indi = [cst_ind[0][i] if i == 0 else cst_ind[0][i-1] + 1 for i in mode_ii]
mode_inde = [cst_ind[0][i] for i in mode_ii]

# calculating fc and Q for the inversion with all events
cst = np.loadtxt(Q_file, skiprows=1)    
print ('allevents\t', cst[mode_indi[0]],'\t', cst[mode_indi[1]])
f_all = ( f0 + (4*np.pi)**-0.5 * cst[mode_indi[0]] )
Q_all = ( 0.5*f_all / (0.5*f0/Q0 + (4*np.pi)**-0.5 * cst[mode_indi[1]]) ) # 12, ordering of coed with cross, cst, dst
fm  = []; lat = []; lon = []; fc = []; Q = []

# calculating the focal mechanism, fc and Q for each event
for i in range(len(cat)):
    # cmt figure
    cmt = cat[i][0].focal_mechanisms[0].moment_tensor.tensor
    fm.append([cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp])
    lat.append(cat[i][0].origins[0].latitude) 
    lon.append(cat[i][0].origins[0].longitude)
    
    # Q figure
    Q_file = "%s/%s/d%s/Iter_10/cst_solver/new_model.damp_%.2e.cst" %(cwd, cmt_id[i], damp, damp)
    cst = np.loadtxt(Q_file, skiprows=1)    
    fc.append( f0 + (4*np.pi)**-0.5 * cst[mode_indi[0]] )
    Q.append( 0.5*fc[0] / (0.5*f0/Q0 + (4*np.pi)**-0.5 * cst[mode_indi[1]]) ) 
    print (cmt_id[i], '\t', cst[mode_indi[0]], '\t', cst[mode_indi[1]])#, fc[i], Q[i]

# Plotting
# CMT figure    
fig = plt.figure()
ax = fig.add_subplot(131)
m = Basemap(projection='moll',lon_0=180,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='lightgray',lake_color='white')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='white')
 
for i in range(len(cat)):
    axins = inset_axes(ax, .1, .1)
    xmax = ax.get_xlim()[1]
    ymax = ax.get_ylim()[1]
    x, y = m(lon[i], lat[i]); h = .17/2; w = .075/2
    x = x/xmax - w/2.; y = y/ymax - h/2.
    ip = InsetPosition(ax, [x, y, w, h])
    axins.set_axes_locator(ip)
    beachball(fm[i], axis=axins)

# Q Figure
#blue: high attenuation 
#red:low attenuation 
ax = fig.add_subplot(132)
m = Basemap(projection='moll',lon_0=180,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='lightgray',lake_color='white')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='white')
dq = Q_all - Q

for i in range(len(cat)):
    x, y = m(lon[i], lat[i])
    #print cmt_id[i], cat[i][0].origins[0].time, q, Q[i], color[i], lat, lon
    scatter = m.scatter(x, y, c=dq[i], vmin=-max(np.abs(dq))*0.75, vmax=max(np.abs(dq))*0.75,
                        s=500, marker="o", alpha=0.9, cmap=cm.seismic, zorder=10)

divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
fig.add_axes(cax)
clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
clb.ax.set_title('dq', position=(0.5, -2))
ax.set_title('$Q_{allevents}=%d$'%(Q_all), fontsize=18)

# fc Figure
ax = fig.add_subplot(133)
m = Basemap(projection='moll',lon_0=180,resolution='c')    
m.drawcoastlines()
m.fillcontinents(color='lightgray',lake_color='white')
m.drawparallels(np.arange(-90.,91.,30.),labels=[True,True,False,False],dashes=[2,2])
m.drawmeridians(np.arange(-180.,181.,60.),labels=[False,False,False,False],dashes=[2,2])
m.drawmapboundary(fill_color='white')
df = f_all - fc

for i in range(len(cat)):
    x, y = m(lon[i], lat[i])
    scatter = m.scatter(x, y, c=df[i], vmin=-max(np.abs(df))*0.9, vmax=max(np.abs(df))*0.9,
                        s=500, marker="o", alpha=0.9, cmap=cm.get_cmap('seismic_r'), zorder=10)

divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="8%", pad=0.1, pack_start=True)
fig.add_axes(cax)
clb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
clb.ax.set_title('df ($\mu Hz$)', position=(0.5, -2))
ax.set_title('$fc_{allevents}=%.2f$'%(f_all), fontsize=18)

plt.suptitle('$%s$'%(mode[mode_name]['name']), fontsize=25)
plt.subplots_adjust(wspace=0, hspace=0, top=0.85, bottom=-0.85, right=0.85, left=-0.85)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])