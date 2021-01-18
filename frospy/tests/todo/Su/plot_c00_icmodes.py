#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:50:20 2019

@author: talavera
"""

import numpy as np
import matplotlib.pyplot as plt
from frospy.core.modes import read as read_modes
from frospy.core.splittingfunc.plot import sens_kernel
from frospy.core.splittingfunc.read import read_cst_AD, read_cst, read_cst_S20RTS
from frospy.util.read import get_mode_names,get_mode_deg
from frospy.core.modes import Modes
from frospy.util.base import split_digit_nondigit
import matplotlib as mpl

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
mpl.rcParams['axes.linewidth'] = 1.8 #set the value globally
chome = '/net/home/talavera/eejit/splitting/'
                         
allmodes = read_modes()
icvar="R"
family = "vs"

if icvar == "vs":
    cwd='%s/04s00/synthetics/IC/QCC/otheric/ic_c00_1Dvs_%smodes.dat' % (chome,
                                                                        family)
    var="v_s"
    plot_num = "b"
    title="shear velocity"
elif icvar == "R":
    cwd = '%s/04s00/synthetics/IC/QCC/otheric/ic_c00_1DR_%smodes.dat' % (chome,
                                                                        family)
    var="R"
    plot_num ="a"
    title="radius"
    
cfile = open(cwd).readlines()
modes = []
modenames = []
modessize = int(len(cfile)/4)
c = 0
k = 0

greens = plt.cm.get_cmap("Greens", 10)
blues = plt.cm.get_cmap("Blues", 10)
#color = [
##         greens(2), greens(5), greens(7),
#         blues(7), blues(5), blues(2), 
#         greens(2), greens(5), greens(7),
#         ]

color = [
         blues(7), blues(5), blues(3), blues(1),
         greens(1), greens(3), greens(5), greens(7),
         ]

mineos = [
          ["2S3",  1259.80],
          ["3S2",  1129.46],
          ["5S2",  2101.46],
          ["6S3",  2844.57],
          ["8S5",  4184.06],
          ["9S3",  3577.73],
          ["9S4",  3899.56],
          ["10S2", 4057.79],
          ["11S2", 4093.68],
          ["11S4", 4782.79],
          ["11S5", 5088.44],
          ["11S6", 5362.27],
          ["13S6", 6173.30],
          ["14S4", 5558.96],
          ["16S5", 6847.56],
          ["16S6", 7165.19],
          ["16S7", 7485.40],
          ["17S8", 7819.09],
          ["18S6", 7973.16],
          ["21S7", 9184.12],
        ]       
mineos = iter(mineos)
fig1 = plt.figure(figsize=(7.48,4.52))
ax = fig1.add_subplot(1,1,1)

for line in cfile:
    line = line.split()
    n = int(line[0])
    mtype = line[1]
    l = int(line[2])
    vs = float(line[3])
    c00 = float(line[4])
    
    if c%4 == 0:
         k = k + 1

    if vs==0.0025:
        m = '%s%s%s' % (n, mtype.upper(), l)
        _m = allmodes.select(m)[0]
        modes.append(m)
        modename = '${}_{%s}%s_{%s}$' % (n, mtype.upper(), l)
        modenames.append(modename)        
        modes_in = ['1', '%s %s %s %s %s' % (n, mtype.upper(), l, l*2, 0)]
        cst_s20, dst = read_cst_S20RTS(modes_in, None)
        c00_s20 = cst_s20[m]['0'][0]

        fc_h = c00+c00_s20
        fc_l = -c00+c00_s20
        fc_s20 = c00_s20
                
#        fc_h = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (c00+c00_s20))
#        fc_l = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (-c00+c00_s20))
#        fc_s20 = ( _m.freq * 1e3 
#                  -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * c00_s20)
#        
        if k == modessize:            
#            ax.scatter(k, fc_s20, label='S20RTS', linewidths=2, color='w', 
#                       edgecolors='r', marker='D', s=60, zorder=25) 
            ax.scatter(k, fc_h, label="$+0.25\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[3], marker='o',s=20)
            ax.scatter(k, fc_l, label="$-0.25\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[4], marker='o',s=20)
        else:
#            ax.scatter(k, fc_s20, linewidths=2, color='w', 
#                       edgecolors='r', marker='D', s=60, zorder=25) 
            ax.scatter(k, fc_h, linewidths=1, color=color[3], marker='o',s=20)
            ax.scatter(k, fc_l, linewidths=1, color=color[4], marker='o',s=20)
       
        if (m != '7S2' and m != '8S2' and m != '10S2' and m != '11S2' and
            m != '1S0' and m != '2S0' and m != '3S0' and m != '4S0'):
            _m = allmodes.select(m)[0]
            cst, dst, cst_e, dst_e = read_cst_AD(modes_in, None)
            c00_AD = cst[m]['0']
            c00_err_AD = cst_e[m]['0']['uncertainty'][0]

            fc_AD = c00_AD
            fc_err_AD = c00_err_AD
            
#            fc_AD = ( _m.freq * 1e3 
#                     -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * c00_AD)
#            fc_err_AD = 1./np.sqrt(4. * np.pi) * c00_err_AD
            
            if k == modessize: 
                ax.errorbar(k, fc_AD, yerr=fc_err_AD, label='Deuss et al., 2013',
                            elinewidth=1, color='k', mfc='w', mec='k', 
                            marker='D', markersize=4,  mew=1, 
                            capsize=2.5, zorder=50)
            else:
               ax.errorbar(k, fc_AD, yerr=fc_err_AD, 
                           elinewidth=1, color='k', mfc='w', mec='k', 
                           marker='D', markersize=4,  mew=1, 
                           capsize=2.5, zorder=50)
        if m == '10S2':
##            GLW
#            fc = ( _m.freq * 1e3 
#                  -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * -14.8464)
#            fc_err = 1./np.sqrt(4. * np.pi) * 0.3152
#            
#            ax.errorbar(k, fc, yerr=fc_err, 
#                        elinewidth=1, color='gray', mfc='w', mec='gray', 
#                        marker='D', markersize=4,  mew=1, 
#                        capsize=2.5, zorder=50, label='GLW')

##            My Inversion            
#            fc = ( _m.freq * 1e3 - 4.0131 * 1e3)
##                  -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * 40.017952)
#            fc_err = 1./np.sqrt(4. * np.pi) * 0.

            fc = 40.017952
            #fc_err = np.array([0.694, 3.376]).reshape((2, 1))
            fc_err = np.array([10, 10]).reshape((2, 1))
            
            ax.errorbar(k, fc, yerr=fc_err, 
                        elinewidth=1, color='r', mfc='w', mec='r', 
                        marker='D', markersize=4,  mew=1, 
                        capsize=2.5, zorder=50, label='Our inversion')

        if m == '11S2':
##            My Inversion 
#            fc = ( _m.freq * 1e3 - 4.0421 * 1e3)
##                  -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * 39.0687065)
#            fc_err = 1./np.sqrt(4. * np.pi) * 0.
 
            fc = 39.0687065
            #fc_err = np.array([7.25e-05, 1.431]).reshape((2, 1))
            fc_err = np.array([10, 10]).reshape((2, 1))
            
            ax.errorbar(k, fc, yerr=fc_err, 
                        elinewidth=1, color='r', mfc='w', mec='r', 
                        marker='D', markersize=4, mew=1, 
                        capsize=2.5, zorder=50, )
        if m == '7S2':
#            fc = ( _m.freq * 1e3 #- 2.5179 * 1e3)
#                  -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * -0.540576041)
#            fc_err = 1./np.sqrt(4. * np.pi) * 0
            fc = -0.540576041
            fc_err = np.array([0.640540299, 0.405428514]).reshape((2, 1))
            
            ax.errorbar(k, fc, yerr=fc_err, 
                        elinewidth=1, color='gray', mfc='w', mec='r', 
                        marker='D', markersize=4,  mew=1, 
                        capsize=2.5, zorder=50,) #label='Talavera-Soza and Deuss, in prep')
 
    elif vs==0.005:
#        fc_h = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (c00+c00_s20))
#        fc_l = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (-c00+c00_s20))

        fc_h = c00+c00_s20
        fc_l = -c00+c00_s20
        if k == modessize: 
            ax.scatter(k, fc_h, 
                       label="$+0.50\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[2], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       label="$-0.50\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[5], marker='o',s=20)
        else:
            ax.scatter(k, fc_h, 
                       linewidths=1, color=color[2], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       linewidths=1, color=color[5], marker='o',s=20)
            
    elif vs==0.0075:
#        fc_h = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (c00+c00_s20))
#        fc_l = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (-c00+c00_s20))

        fc_h = c00+c00_s20
        fc_l = -c00+c00_s20
        if k == modessize: 
            ax.scatter(k, fc_h, 
                       label="$+0.75\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[1], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       label="$-0.75\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[6], marker='o',s=20)
        else:
            ax.scatter(k, fc_h, 
                       linewidths=1, color=color[1], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       linewidths=1, color=color[6], marker='o',s=20)                      
    elif vs==0.01:
#        f_mineos = mineos.next()[1]
#        c00_mineos = (f_mineos - _m.freq * 1e3) * np.sqrt(4. * np.pi)
        
#        fc_h = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (c00+c00_s20))
#        fc_l = ( _m.freq * 1e3 
#              -_m.freq * 1e3 + 1./np.sqrt(4. * np.pi) * (-c00+c00_s20))
        
        fc_h = c00+c00_s20
        fc_l = -c00+c00_s20
#        print m, c00+c00_s20, c00_mineos + c00_s20, c00_s20, _m.freq * 1e3
#        c00_mineos = c00_mineos + c00_s20
        if k == modessize: 
            ax.scatter(k, fc_h, 
                       label="$+1.00\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[0], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       label="$-1.00\%$ ${{{0}}}$".format(var),# with S20RTS",
                       linewidths=1, color=color[7], marker='o',s=20)
        else:
#            ax.scatter(k, c00_mineos, 
#                       linewidths=1, color="b", marker='s',s=20)
            ax.scatter(k, fc_h, 
                       linewidths=1, color=color[0], marker='o',s=20)
            ax.scatter(k, fc_l, 
                       linewidths=1, color=color[7], marker='o',s=20)
#        print(m, k, fc_l,-c00,c00_s20)
    c = c + 1

ax.set_xticks(np.linspace(1, modessize, modessize))
ax.set_xticklabels(modenames, rotation='horizontal', fontsize=9)
ax.set_xlim(0.5,k+0.5)
#ax.set_ylabel("$f_c$ ($\mu$Hz)", fontsize=10)
ax.set_ylabel("$c_{00}$ ($\mu$Hz)", fontsize=10)
ax.tick_params(axis='y', which='both', labelsize=6)
#plt.yticks(fontsize=8)
plt.axhline(y=0, color='lightgray', 
            linestyle='-',  linewidth=1.5, zorder=0)

# legend
handles, labels = ax.get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
elabels, ehandles = zip(*sorted(zip(labels[0:4], handles[0:4]), 
                                key=lambda t: t[0], reverse=True))
labels, handles = elabels+labels[4::], ehandles+handles[4::]
ehandles = [h[0] for h in handles[8::]]
for h in ehandles: h.set_linestyle("")
handles = handles[0:8] + tuple(ehandles)
ax.legend(handles, labels, frameon=False, 
          loc="lower left", fontsize=8,
          bbox_to_anchor=(-0.01, -0.018),
          columnspacing=0,handletextpad=0.0001,
          labelspacing=0.25)

ax.set_title("{0}) Varying 1D inner core {1} ${{{2}}}$".format(plot_num,
                                                                title,var))
#ax.set_ylim([-104.69499346927807, 99.8741231092651])
#fig1
fig1.savefig('c00_ic%smodes_%schange'%(family, icvar),dpi=400,
              bbox_inches='tight',
              pad_inches=0.05,)

#bbox = ax.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
#x, y = bbox.width, bbox.height

####
#
#mpl.rcParams['axes.linewidth'] = 1. #set the value globally
#
#width = np.ones(modessize)*1.25
#gridspec_kw = {'width_ratios': width, 'wspace':0}
#fig2, axes = plt.subplots(nrows=1, ncols=modessize, sharey=True, sharex=True,
#                         gridspec_kw=gridspec_kw, )
#
#fig2.set_size_inches(7.48, 1.13)
#legend_show = False
#text= True
#for ax, mode, name in zip(axes.flat, modes, modenames):   
#    m = allmodes.select(mode)
#    sens_kernel(m[0], title=None, ax=ax, legend_show=legend_show)
#    if text:
#        ax.text(-6.5,5700,'TZ', fontsize=5.)
#        ax.text(-6.5,3400,'CMB',fontsize=5.)
#        ax.text(-6.5,1220,'ICB',fontsize=5.)
#        text=False
#        ax.spines['left'].set_linewidth(1.8)
#    for axis in ['top','bottom']:
#        ax.spines[axis].set_linewidth(1.8)
#    legend_show = False    
#ax.spines['right'].set_linewidth(1.8)
##ax.set_ylim(6000,0)
##ax.invert_yaxis()
#ax.axes.xaxis.set_ticklabels([])
#ax.axes.yaxis.set_ticklabels([])
#ax.axes.get_xaxis().set_ticks([])
#ax.axes.get_yaxis().set_ticks([])
#fig2
#fig2.savefig('icmodes_sens_3',dpi=400,
#              bbox_inches='tight',
#              pad_inches=0.,)

#name = '7s2'
#m = read_modes(modenames=name)[0]
#sens_kernel(m, title=name)