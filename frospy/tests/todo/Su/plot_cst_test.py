#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 16:30:27 2018

@author: talavera
"""
import numpy as np
from obspy.core.util.attribdict import AttribDict
from frospy.postprocessing.AnalyseRuns.plot import misfits_cst
from frospy.util.read import read_cst, read_modes_in, get_cst_dst
from frospy.util.base import (sc_degrees, cc_degrees)
from frospy.core.splittingfunc import get_covariances
from frospy.util.read import read_cst, read_cst_dst, read_modes_in, get_cst_dst
from frospy.util.nmfile import sort_AttribDict
from frospy.util.base import (sc_degrees, cc_degrees)
import matplotlib as mpl
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import (FormatStrFormatter,
                               MaxNLocator)
from frospy.core.splittingfunc import plot_cst, get_cst_list, get_cst_d00_list
import os

#folder='/net/home/talavera/radial-inversion/05s00/05s00/deepevents/cst+d00/prem'
#folder='/net/home/talavera/radial-inversion/02s00/02s00-07s02/deepevents/cst+d00/prem'
folder = os.getcwd()
modes_dir='%s/d0.001/config'%folder
cst_file='%s/d0.001/Iter_10/cst_solver/mcst.dat'%folder


#misfits_cst(folder, "Z", "damping", startiter=1, enditer=2, plot_coeff=True)
modes, modes_cc = read_modes_in(modes_dir)

#with open(cst_file, 'r') as fh:
#    c_tmp = np.genfromtxt(fh)
#    if c_tmp.ndim > 1:
#        # noc is set to False, if Arwens format is read
#        noc = False
#        c = c_tmp.transpose()[0]
#    else:
#        # noc is set to the number of cst, if Haydars format is read
#        noc = int(c_tmp[0])
#        c = c_tmp[1:]

def get_dst_list(infile, modes_dir, degree):

    dst_list = []
    for fh, md in zip(infile, modes_dir):
        cst, dst = read_cst_dst(fh, md)
        dst_sorted = AttribDict()
        for k in dst.keys():
            dst_tmp = sort_AttribDict(dst[k], 'int')
            dst_sorted[k] = dst_tmp
        dst_list.append(dst_sorted)

    if degree is not 'all':
        for dst in dst_list:
            for name, coeff in dst.iteritems():
                new_coeff = []
                for deg in coeff:
                    if int(deg.keys()[0]) in degree:
                        new_coeff.append(deg)
                dst[name] = new_coeff
    return dst_list
cst_deg='all'
dst_list = get_dst_list([cst_file], [modes_dir], cst_deg)

def get_cst_d00_list(infile, modes_dir, degree):
    cst_list = []
    d00_list = []
    models = ['RR', 'TZ', 'S20RTS', 'REM']
    for fh, md in zip(infile, modes_dir):
        cst, dst = read_cst(fh, md)

        cst_sorted = AttribDict()
        for k in cst.keys():
            cst_tmp = sort_AttribDict(cst[k], 'int')
            cst_sorted[k] = cst_tmp
        cst_list.append(cst_sorted)

        d00_sorted = AttribDict()
        for k in dst.keys():
            if k.find("-") == -1 and fh not in models: # only d0 of self modes
                d00_tmp = sort_AttribDict(dst[k], 'int')
                d00_sorted[k] = d00_tmp
        d00_list.append(d00_sorted)

    if degree is not 'all':
        for cst in cst_list:
            for name, coeff in cst.iteritems():
                new_coeff = []
                for deg in coeff:
                    if int(deg.keys()[0]) in degree:
                        new_coeff.append(deg)
                cst[name] = new_coeff

        for dst in d00_list:
            for name, coeff in dst.iteritems():
                new_coeff = []
                for deg in coeff:
                    if int(deg.keys()[0]) in degree:
                        new_coeff.append(deg)
                dst[name] = new_coeff
    return cst_list, d00_list

cst_deg='all'
#cst_list, d00_list = get_cst_d00_list([cst_file], [modes_dir], cst_deg)

cst_files=[]
modes_dir=[]
labellist=[]
#damp=['d0.00001', 'd0.0001', 'd0.001', 'd0.01', 'd0.1', 'd1', 'd10', 'd100', 'd1000', 'd10000']
damp=['d0.00001','d10000']
for d in damp:
    path = "%s/%s/Iter_10/cst_solver/mcst.dat" % (folder, d)
    cst_files.append(path)
    modes_dir.append('%s/%s/config'%(folder,d))
    labellist.append(d)
    
cst_list, d00_list = get_cst_d00_list(cst_files, modes_dir, cst_deg)

plot_cst(modes_dir, cst_files, covmatrix=None, labellist=labellist,
         plot=True, cst_deg='all')


def plot_dst(modes_dir, dst_file=None, covmatrix=None, labellist=None,
             dst_deg='all', iter_no='10',
             save_fig=False, yscale='auto', plot=True, verbose=False):
    """
    param folder: path/s to inversion files

    type folder: string

    param modes_dir: path to directory containing:folder
                         modes.in
                         modes_cc.in
                     as defined for synseis hybrid

    type  modes_dir: string

    param degree: list of degrees to be plotted
    type  degree: list of integers or 'all'

    param yscale: 'auto' or 'global'
    type  yscale: string

    param cstfile: file containing cst coefficients

    Example for new coefficients of hybrid code:
    Modes 0S15 - 0T16
    The order of coefficients are as follows:
    indices    mode
    1-231      0 S 15 cst
    232-237    0 T 16 cst
    238-242    0 S 15 - 0 T 16 cst
    243        0 S 15 d00
    244        0 T 16 d00
    """

    modes = AttribDict()
    
    if save_fig:
        mpl.rcParams.update({'font.size': 15})

    infiles = []
    if type(cst_file) != list:
        dst_file = [dst_file]
    for f in cst_file:
        infiles.append(f)

    if type(modes_dir) != list:
        modes_dir = [modes_dir]

    if dst_deg != 'all':
        for i, c in enumerate(dst_deg):
            dst_deg[i] = int(c)

    model_count = 0
    if modes_dir is not None:
        for model in ['RR', 'TZ', 'S20RTS', 'REM']:
            if model in infiles:
                modes_dir.append(modes_dir[0])
                model_count += 1
        dst_list = get_dst_list(infiles, modes_dir, dst_deg)
        if covmatrix is not None:
            cst_cov, dst_cov = get_covariances(covmatrix, modes_dir)
            
    # getting dst only
    d = [AttribDict() for i in range(len(dst_list))]
    for _i, dst in enumerate(dst_list):
        for name, coeffs in dst.iteritems():
          c=[]
          for cval in coeffs:
                degree = cval.keys()[0]
                if not degree=='0': 
                    c.append(cval)
                    d[_i].update({name:c})
    dst_list = d                
    if model_count != 0:
        colormap = cm.rainbow(np.linspace(0, 1, len(dst_list)-model_count))
        colormap = np.vstack((colormap, 0.5 * np.ones(colormap[-1].shape)))
    else:
        colormap = cm.rainbow(np.linspace(0, 1, len(dst_list)))

    for _i, dst in enumerate(dst_list):
        if verbose:
            print("Plotting coefficients for %s" % infiles[_i])
        for name, coeffs in dst.iteritems():
            a_i = 0
            N = len(coeffs)
            if name not in modes:
                modes[name] = AttribDict()
                fig = plt.figure()
                fig.suptitle(name)
                modes[name]['fig'] = fig
                modes[name]['axes'] = AttribDict()

            else:
                fig = modes[name]['fig']

            axes = modes[name]['axes']
            # init plots
            for cval in coeffs:
                degree = cval.keys()[0]

                if degree not in axes:
                    axes[degree] = AttribDict()
                    ylim = 0
                    xlim = 0
                    axes[degree]['xmax'] = 0
                    axes[degree]['ymax'] = 0
                    axes[degree]['ax'] = plt.subplot2grid((N, 2), (a_i, 0),
                                                          colspan=2)

                    axes[degree]['ax'].set_title(r"$s=%i$" % int(degree))
                    a_i += 1

                if abs(cval.values()[0]).max() > abs(ylim):
                    ylim = abs(cval.values()[0]).max()
                    axes[degree]['ymax'] = ylim
                if len(cval.values()[0]) > xlim:
                    xlim = len(cval.values()[0])
                    axes[degree]['xmax'] = xlim

            for _j, cval in enumerate(coeffs):
                degree = cval.keys()[0]
                y = cval.values()[0]
                x = np.arange(len(y))
                ax = axes[degree].ax
                xlim = axes[degree].xmax
                ylim = axes[degree].ymax
                ax.hlines(0, -0.5, xlim+1, linestyle='-', color='black',
                          linewidth=.45)
                if labellist is not None:
                    cstlabel = labellist[_i]
                    if cstlabel == 'RR':
                        marker = 'o'
                        color = 'grey'
                        linestyle = ':'
                    elif cstlabel == 'REM':
                        marker = 'o'
                        color = 'grey'
                        linestyle = '--'
                    elif cstlabel == 'TZ':
                        marker = 'o'
                        color = 'grey'
                        linestyle = '-.'
                    else:
                        color = colormap[_i]
                        marker = '^'
                        linestyle = '-'
                    if covmatrix:
                        try:
                            yerr = cst_cov[_i][name][_j].values()[0]
                            ax.errobar(x, y, yerr=yerr, marker='^',
                                       label=infiles[_i].split('/')[-1],
                                       linewidth=0.8, barsabove=True,
                                       elinewidth=0.5,
                                       ecolor=color,
                                       capsize=5,
                                       capthick=0.5,
                                       color=color)
                        except Exception:
                            ax.plot(x, y, marker=marker, label=cstlabel,
                                    linewidth=0.8, linestyle=linestyle,
                                    color=color)
                    else:
                        ax.plot(x, y, marker=marker, label=cstlabel,
                                linewidth=0.8, linestyle=linestyle,
                                color=color)
                else:
                    if covmatrix:
                        try:
                            yerr = cst_cov[_i][name][_j].values()[0]

                            ax.errorbar(x, y, yerr=yerr, marker='^',
                                        label=infiles[_i].split('/')[-1],
                                        linewidth=0.8, barsabove=True,
                                        elinewidth=0.5,
                                        ecolor=colormap[_i],
                                        capsize=5,
                                        capthick=0.5,
                                        color=colormap[_i])
                        except Exception:
                            ax.plot(x, y, marker='^',
                                    label=infiles[_i].split('/')[-1],
                                    linewidth=0.8,
                                    color=colormap[_i])
                    else:
                        ax.plot(x, y, marker='^',
                                label=infiles[_i].split('/')[-1],
                                linewidth=0.8,
                                color=colormap[_i])

                ax.xaxis.set_major_locator(MaxNLocator(integer=True))

                ax.set_xlim((-0.5, xlim+1))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    legend_set = False
    for ax in axes.itervalues():
        if len(ax.ax.lines) == len(dst_list):
            ax.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            legend_set = True
            break
        if not legend_set:
            if cst_deg == 'all':
                axes['0'].ax.legend(loc='center left',
                                    bbox_to_anchor=(1, 0.5))
            else:
                axes[str(cst_deg[0])].ax.legend(loc='center left',
                                                bbox_to_anchor=(1, 0.5))

    if save_fig:
        for fi in plt.get_fignums():
            plt.figure(fi)
            filename = str(fi) + '_dst.png'
            fig = plt.gcf()
            fig.set_size_inches(12, 8)
            plt.savefig(filename, dpi=300, orientation='landscape')
            print("Saving figure %s" % filename)
        plt.close('all')
    elif plot:
        plt.show()
    return

def plot_cst(modes_dir, cst_file=None, covmatrix=None, labellist=None,
             cst_deg='all', iter_no='10',
             save_fig=False, yscale='auto', plot=True, verbose=False):
    """
    param folder: path/s to inversion files

    type folder: string

    param modes_dir: path to directory containing:folder
                         modes.in
                         modes_cc.in
                     as defined for synseis hybrid

    type  modes_dir: string

    param degree: list of degrees to be plotted
    type  degree: list of integers or 'all'

    param yscale: 'auto' or 'global'
    type  yscale: string

    param cstfile: file containing cst coefficients

    Example for new coefficients of hybrid code:
    Modes 0S15 - 0T16
    The order of coefficients are as follows:
    indices    mode
    1-231      0 S 15 cst
    232-237    0 T 16 cst
    238-242    0 S 15 - 0 T 16 cst
    243        0 S 15 d00
    244        0 T 16 d00
    """
    cst_file=['%s/d10000/Iter_10/cst_solver/mcst.dat'%folder,
              '%s/d0.00001/Iter_10/cst_solver/mcst.dat'%folder]

    modes_dir=['%s/d10/config'%folder,
               '%s/d0.001/config'%folder
               ]
#    modes_dir=modes_dir[0]
    
    covmatrix=None
    labellist=None
    cst_deg='all'
    iter_no='10'
    save_fig=False
    yscale='auto'
    plot=True
    verbose=False
    modes = AttribDict()

    if save_fig:
        mpl.rcParams.update({'font.size': 15})

    infiles = []
    if type(cst_file) != list:
        cst_file = [cst_file]
    for f in cst_file:
        infiles.append(f)

    if type(modes_dir) != list:
        modes_dir = [modes_dir]

    if cst_deg != 'all':
        for i, c in enumerate(cst_deg):
            cst_deg[i] = int(c)

    model_count = 0
    if modes_dir is not None:
        for model in ['RR', 'TZ', 'S20RTS', 'REM']:
            if model in infiles:
                modes_dir.append(modes_dir[0])
                model_count += 1
        # cst_list = get_cst_list(infiles, modes_dir, cst_deg)
        cst_list, d00_list = get_cst_d00_list(infiles, modes_dir, cst_deg)
        if covmatrix is not None:
            cst_cov, dst_cov = get_covariances(covmatrix, modes_dir)
            
    d = [AttribDict() for i in range(len(dst_list))]
    for _i, dst in enumerate(dst_list):
        for name, coeffs in dst.iteritems():
          for cval in coeffs:
                degree = cval.keys()[0]
                if not degree=='0':   
                    d[_i].update({name:cval})
                
    if model_count != 0:
        colormap = cm.rainbow(np.linspace(0, 1, len(cst_list)-model_count))
        colormap = np.vstack((colormap, 0.5 * np.ones(colormap[-1].shape)))
    else:
        colormap = cm.rainbow(np.linspace(0, 1, len(cst_list)))

    for _i, cst in enumerate(cst_list):
        print("Plotting coefficients for %s" % infiles[_i])
        for name, coeffs in cst.iteritems():
            a_i = 0
            N = len(coeffs)
            if name not in modes:
                modes[name] = AttribDict()
                fig = plt.figure()
                fig.suptitle(name)
                modes[name]['fig'] = fig
                modes[name]['axes'] = AttribDict()

            else:
                fig = modes[name]['fig']

            axes = modes[name]['axes']
            # init plots
            for cval in coeffs:
                degree = cval.keys()[0]

                if degree not in axes:
                    axes[degree] = AttribDict()
                    ylim = 0
                    xlim = 0
                    axes[degree]['xmax'] = 0
                    axes[degree]['ymax'] = 0
                    axes[degree]['ax'] = plt.subplot2grid((N, 2), (a_i, 0),
                                                          colspan=2)

                    axes[degree]['ax'].set_title(r"$s=%i$" % int(degree))
                    # trans = axes[degree]['ax'].get_xaxis_transform()
                    # axes[degree]['ax'].annotate(r"$s=%i$" % int(degree),
                    #                             xy=(xlim/2., .84),
                    #                             xycoords=trans)
                    a_i += 1

                if abs(cval.values()[0]).max() > abs(ylim):
                    ylim = abs(cval.values()[0]).max()
                    axes[degree]['ymax'] = ylim
                if len(cval.values()[0]) > xlim:
                    xlim = len(cval.values()[0])
                    axes[degree]['xmax'] = xlim

            for _j, cval in enumerate(coeffs):
                degree = cval.keys()[0]
                y = cval.values()[0]
                if degree=='0':
                    for name_d00, coeffs_d00 in d00_list[_i].iteritems():
                        if name==name_d00:
                            y = np.append(y, coeffs_d00[0]['0'])
                            print(name_d00,coeffs_d00[0]['0'], y) 

                x = np.arange(len(y))
                #print(name, degree, x, y)
                ax = axes[degree].ax
                xlim = axes[degree].xmax
                ylim = axes[degree].ymax
                ax.hlines(0, -0.5, xlim+1, linestyle='-', color='black',
                          linewidth=.45)
                if labellist is not None:
                    cstlabel = labellist[_i]
                    if cstlabel == 'RR':
                        marker = 'o'
                        color = 'grey'
                        linestyle = ':'
                    elif cstlabel == 'REM':
                        marker = 'o'
                        color = 'grey'
                        linestyle = '--'
                    elif cstlabel == 'TZ':
                        marker = 'o'
                        color = 'grey'
                        linestyle = '-.'
                    else:
                        color = colormap[_i]
                        marker = '^'
                        linestyle = '-'
                    if covmatrix:
                        try:
                            yerr = cst_cov[_i][name][_j].values()[0]
                            ax.errobar(x, y, yerr=yerr, marker='^',
                                       label=infiles[_i].split('/')[-1],
                                       linewidth=0.8, barsabove=True,
                                       elinewidth=0.5,
                                       ecolor=color,
                                       capsize=5,
                                       capthick=0.5,
                                       color=color)
                        except Exception:
                            ax.plot(x, y, marker=marker, label=cstlabel,
                                    linewidth=0.8, linestyle=linestyle,
                                    color=color)
                    else:
                        ax.plot(x, y, marker=marker, label=cstlabel,
                                linewidth=0.8, linestyle=linestyle,
                                color=color)
                else:
                    if covmatrix:
                        try:
                            yerr = cst_cov[_i][name][_j].values()[0]

                            ax.errorbar(x, y, yerr=yerr, marker='^',
                                        label=infiles[_i].split('/')[-1],
                                        linewidth=0.8, barsabove=True,
                                        elinewidth=0.5,
                                        ecolor=colormap[_i],
                                        capsize=5,
                                        capthick=0.5,
                                        color=colormap[_i])
                        except Exception:
                            ax.plot(x, y, marker='^',
                                    label=infiles[_i].split('/')[-1],
                                    linewidth=0.8,
                                    color=colormap[_i])
                    else:
                        ax.plot(x, y, marker='^',
                                label=infiles[_i].split('/')[-1],
                                linewidth=0.8,
                                color=colormap[_i])

                ax.xaxis.set_major_locator(MaxNLocator(integer=True))

                ax.set_xlim((-0.5, xlim+1))

                # if degree != '0':
                #     if yscale == 'global':
                #         ax.set_ylim((-1.15 * ylim, 1.15 * ylim))
                #
                #     elif yscale == 'auto':
                #         limit = list(ax.get_ylim())
                #         if y.min() <= .85 * limit[0] or y.min() == 0:
                #             if y.min() == 0.:
                #                 limit[0] = -1.15
                #             else:
                #                 limit[0] = 1.15 * y.min()
                #
                #         if y.max() >= .85 * limit[1] or y.max() == 0:
                #             if y.max() == 0.:
                #                 limit[1] = 1.15
                #             else:
                #                 limit[1] = 1.15 * y.max()
                #         ax.set_ylim(tuple(limit))

                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    legend_set = False
    for ax in axes.itervalues():
        if len(ax.ax.lines) == len(cst_list):
            ax.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            legend_set = True
            break
        if not legend_set:
            if cst_deg == 'all':
                axes['0'].ax.legend(loc='center left',
                    bbox_to_anchor=(1, 0.5))
            else:
                axes[str(cst_deg[0])].ax.legend(loc='center left',
                     bbox_to_anchor=(1, 0.5))

    # max_deg = max(axes.keys())
    # for deg, values in axes.iteritems():
    #     if deg != max_deg:
    #         values.ax.get_xaxis().set_ticks([])

    if save_fig:
        for fi in plt.get_fignums():
            plt.figure(fi)
            filename = str(fi) + '_cst.png'
            fig = plt.gcf()
            fig.set_size_inches(12, 8)
            plt.savefig(filename, dpi=300, orientation='landscape')
            print("Saving figure %s" % filename)
        plt.close('all')
    elif plot:
        plt.show()
    return