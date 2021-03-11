# -*- coding: utf-8 -*-
"""
Module for plotting nmPy SplittingFunc objects .

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import, print_function
from obspy.core.util import AttribDict

import numpy as np
# import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import cm
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import (FormatStrFormatter, MaxNLocator)

from frospy.util.base import split_digit_nondigit, sort_human, update_progress
from frospy.util.base import max_cc_degrees, find_unique_name
from frospy.plot.nmplt import shiftedColorMap, get_iter_colormap
import frospy.data.basemap.plates as tplates
import frospy.data.basemap.llsvp as LLSVP

import os
import glob
import subprocess


class Bin(object):
    def __init__(self):
        self.sc_cstkernels = "/net/home/talavera/bin/mdcplmrho_all_cstkernels"
        self.sc_dstkernels = "/net/home/talavera/bin/mdcplmrho_kmr_cstkernels"
        self.cc_kernels = "/net/home/talavera/bin/mdcplmrho_allC_cstkernels"


def plot_cst_partial(path, labels=None, weights=None):
    fig = Figure()
    ax  = fig.add_subplot(111)

    if labels is None:
        # num = np.arange(1, len(path)+1).astype(str)
        # name = np.chararray(len(path), itemsize=11)
        # name[:] = 'cst '
        # labels = np.core.defchararray.add(name, num)
        labels = path

    if weights is None:
        weights = np.ones(len(path))
    for i, f in enumerate(path):
        y = np.genfromtxt(f)
        y = weights[i] * y.diagonal()
        x = range(y.size)


        ax.plot(x, y, label=labels[i])
    ax.legend()
    # plt.show()

    return


def sens_kernel_branch(modes, colormap='rainbow', kernel='vs', savefig=False,
                       newfigure=False, yticks=True):
    """
    To produce an animation use this command:
    mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts \
    vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
    """
    colormap = get_iter_colormap(modes, colormap)
    cmap = next(colormap)
    fig, ax = sens_kernel(modes[0], color=[cmap, cmap, cmap], kernel=kernel)
    if savefig is True:
        fig.savefig('mode001.png')
    for i, m in enumerate(modes[1:]):
        cmap = next(colormap)
        if newfigure is True:
            ax = fig = None
        fig, ax = sens_kernel(m, ax=ax, fig=fig, color=[cmap, cmap, cmap],
                              kernel=kernel)
        if savefig is True:
            fig.savefig('mode%.3d.png' % int(i+2))
            update_progress((i+1)/float(len(modes[1:])), 'Printing')
    return fig, ax


def sens_kernel(mode, ax=None, fig=None, title=True, show=False, savefig=False,
                legend_show=True, color='auto', kernel='all',
                kind="cst", **kwargs):
    """
    :params mode: frospy.core.modes.Mode
    """
    bins = Bin()
    if not os.path.exists(bins.sc_cstkernels):
        raise OSError('senskernel bin not found!')

    if 'fontsize' in kwargs:
        fontsize = float(kwargs['fontsize'])
    else:
        fontsize = 8

    if 'linewidth' in kwargs:
        linewidth = float(kwargs['linewidth'])
    else:
        linewidth = 1

    if glob.glob('*-kernel.dat'):
        os.system('rm *-kernel.dat')

    if ax is None:
        fig = Figure()
        ax  = fig.add_subplot(111)
        fig.set_size_inches(4, 11)
        ax.lines = []

    ## get the kernel.dat files with modeplotqt
    #os.system('echo "j" > input')
    #os.system('echo %03d %s %03d >> input' % (mode.n,  mode.type,  mode.l))
    #res = subprocess.Popen('/net/home/jagt/modes/kernels/bin/modeplotqt <'
    #                       ' input', shell=True, stdout=subprocess.PIPE,
    #                       stderr=subprocess.PIPE)
    #os.system('rm input')

    #output, error = res.communicate()
    #if res.returncode != 2:
    #    raise Exception(error)
    #    print('error: %s') % res.returncode

    ## get the kernel.dat files with mdcplmrho, works for all degrees
    ## e.g. 0s2 flips the sign of its kernel between deg=0(+) and deg=2(-)
    if kind == "cst":
        # sc_cstkernels = "/net/home/talavera/bin/mdcplmrho_all_cstkernels"
        sc_cstkernels = bins.sc_cstkernels
        ktypes = ["s", "p", "r"]
        vp = []
        vs = []
        rho = []
        depth = []
    elif kind == "dst":
        # sc_cstkernels = "/net/home/talavera/bin/mdcplmrho_kmr_cstkernels"
        sc_cstkernels = bins.sc_dstkernels
        ktypes = ["k", "m"]
        kappa = []
        mu = []
        depth = []

    cwd = os.getcwd()
    tmp_path = 'tmp_sckernel'

    if os.path.exists(tmp_path):
        tmp_path = find_unique_name(tmp_path)
    os.makedirs(tmp_path)
    os.chdir(tmp_path)

    modes = []
    mtype = []
    os.system('echo "1" > modes.in')
    os.system('echo "%03d %s %03d" >> modes.in' % (mode.n,
                                                   mode.type.lower(),
                                                   mode.l))

    if mode.name == '0S2':
        cdeg = 2
    else:
        cdeg = 0
    # kernels generation
    for k in ktypes:
        os.system('echo %s > input' % k)
        os.system('echo "%s" >> input' % cdeg) # calculated for deg=2

        res = subprocess.Popen('%s < input' % sc_cstkernels, shell=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        output, error = res.communicate()

        if k == "s":
            kname = "vs"
        if k == "p":
            kname = "vp"
        if k == "r":
            kname = "rho"
        if k == "m":
            kname = "mu"
        if k == "k":
            kname = "kappa"
        os.system('mv kernel.dat %s-kernel.dat' % kname)
        os.system('rm ast* cst.dat mcst.dat mdcpl.out raw.dat')
    os.system('rm input modes.in')

    if res.returncode != 0:
        raise Exception(error)
        print('error: %s') % res.returncode

    if kind == "cst":
        with open('rho-kernel.dat') as f:
            for line in f:
                cols = line.split()
                #depth.append(float(cols[0])) # modeplotqt
                #rho.append(float(cols[1])) # modeplotqt
                depth.append(float(cols[1])) # mdcplmrho
                rho.append(float(cols[0])) # mdcplmrho
        os.system('rm rho-kernel.dat')
        rhomax = np.amax(np.abs(rho))

        with open('vs-kernel.dat') as f:
            for line in f:
                cols = line.split()
                #vs.append(float(cols[1])) # modeplotqt
                vs.append(float(cols[0])) # mdcplmrho
        os.system('rm vs-kernel.dat')
        vsmax = np.amax(np.abs(vs))
        intmax = np.amax([rhomax, vsmax])

        if mode.type == 'S' and kernel in ('all', 'vp'):
            with open('vp-kernel.dat') as f:
                for line in f:
                    cols = line.split()
                    #vp.append(float(cols[1])) # modeplotqt
                    vp.append(float(cols[0])) # mdcplmrho
            os.system('rm vp-kernel.dat')
            vpmax = np.amax(np.abs(vp))
            intmax = np.amax([intmax, vpmax])
    elif kind == "dst":
        with open('mu-kernel.dat') as f:
            for line in f:
                cols = line.split()
                depth.append(float(cols[1])) # mdcplmrho
                mu.append(float(cols[0])) # mdcplmrho
        os.system('rm mu-kernel.dat')
        mumax = np.amax(np.abs(mu))

        if mode.type == 'S' and kernel in ('all', 'kappa'):
            with open('kappa-kernel.dat') as f:
                for line in f:
                    cols = line.split()
                    kappa.append(float(cols[0])) # mdcplmrho
            os.system('rm kappa-kernel.dat')
            kappamax = np.amax(np.abs(kappa))
            intmax = np.amax([mumax, kappamax])

    ax.set_ylim(0, 6400)
    ax.set_xlim(-1.1*intmax, 1.1*intmax)

    ax.axhline(5700, color='k', lw=0.1)  # TZ
    ax.axhline(3480, color='k', lw=0.1)  # CMB
    ax.axhline(1220, color='k', lw=0.1)  # ICB
    ax.axvline(0, color='k', lw=0.1)

    if color == 'auto':
        rho_clr = 'grey'
        vs_clr = 'red'
        vp_clr = 'k'
        mu_clr = 'red'
        kappa_clr = 'k'
    else:
        rho_clr = color[0]
        vs_clr = color[1]
        vp_clr = color[2]
        mu_clr = color[1]
        kappa_clr = color[2]

    legend_ax = []
    legend = []
    if kind == "cst":
        if kernel in ('all', 'rho'):
            han1, = ax.plot(rho, depth, color=rho_clr, linewidth=linewidth)
            legend_ax.append(han1)
            legend.append(r'$\rho$')

        if kernel in ('all', 'vs'):
            han2, = ax.plot(vs, depth, color=vs_clr, linewidth=linewidth)
            legend_ax.append(han2)
            legend.append(r'$v_s$')

        if mode.type == 'S' and kernel in ('all', 'vp'):
            han3, = ax.plot(vp, depth, color=vp_clr, linewidth=linewidth)
            legend_ax.append(han3)
            legend.append(r'$v_p$')

    elif kind == "dst":
        if kernel in ('all', 'mu'):
            han2, = ax.plot(mu, depth, color=mu_clr, linewidth=linewidth)
            legend_ax.append(han2)
            legend.append(r'$q_\mu$')

        if mode.type == 'S' and kernel in ('all', 'kappa'):
            han3, = ax.plot(kappa, depth, color=kappa_clr, linewidth=linewidth)
            legend_ax.append(han3)
            legend.append(r'$q_\kappa$')

    if 'ticks' in kwargs:
        ticks = kwargs['ticks']
    else:
        ticks = True



    if legend_show and ticks:
        if 'bbox_to_anchor' in kwargs:
            _bbox = kwargs['bbox_to_anchor']
        else:
            _bbox = (-0.08, -0.03)
        ax.legend(legend_ax, legend, frameon=False,
                  bbox_to_anchor=_bbox,
                  handlelength=0.5, handletextpad=0.1, loc='lower left',
                  fontsize=fontsize)

    if legend_show and not ticks:
        if 'bbox_to_anchor' in kwargs:
            _bbox = kwargs['bbox_to_anchor']
        else:
            _bbox = (0.6, 1.2)
        ax.legend(legend_ax, legend, ncol=3,
                  handlelength=0.3, handletextpad=0.1, columnspacing=0.25,
                  loc='upper center', bbox_to_anchor=_bbox, frameon=False,
                  fontsize=fontsize, borderpad=0)


    if not ticks:
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.text(-intmax, 5700, "TZ",  fontsize=fontsize*0.8, va="bottom")
        ax.text(-intmax, 3480, "CMB", fontsize=fontsize*0.8, va="bottom")
        ax.text(-intmax, 1220, "ICB", fontsize=fontsize*0.8, va="bottom")

    if title:
        title = r"$_{%s}%s_{%s}$" % (mode.n,  mode.type,  mode.l)
        ax.set_title(title, fontsize=fontsize*1.5, va="center")
    if show:
        plt.show()
    if savefig and ax is None:
        fig.savefig('%s_kernel.png' % mode.name, orientation='landscape',
                    dpi=400)

    #os.system('rm @imagen') # modeplotqt
    #os.system('rm eig.dat') # modeplotqt
    os.chdir(cwd) # mdcplmrho
    # os.removedirs(tmp_path) # mdcplmrho
    os.system('rm -rf %s' % tmp_path)
    return fig, ax


def sensC_kernel(mode, ax=None, fig=None, title=True, show=False, savefig=False,
                legend_show=True, color='auto', kernel='all', **kwargs):
    """
    :params mode: frospy.core.modes.Mode
    """
    bins = Bin()
    if not os.path.exists(bins.cc_kernels):
        raise OSError('cc sens kernel bin not found!')
    if 'fontsize' in kwargs:
        fontsize = float(kwargs['fontsize'])
    else:
        fontsize = 8

    if 'linewidth' in kwargs:
        linewidth = float(kwargs['linewidth'])
    else:
        linewidth = 1

    if glob.glob('*-kernel.dat'):
        os.system('rm *-kernel.dat')

    if ax is None:
        fig = Figure()
        ax  = fig.add_subplot(111)
        fig.set_size_inches(4, 11)
        ax.lines = []

    cc_cstkernels = bins.cc_kernels
    # cc_cstkernels = "/net/home/talavera/bin/mdcplmrho_allC_cstkernels"

    cwd = os.getcwd()
    tmp_path = 'tmp_ckernel'

    if os.path.exists(tmp_path):
        tmp_path = find_unique_name(tmp_path)
    os.makedirs(tmp_path)
    os.chdir(tmp_path)

    modes = []
    mtype = []
    os.system('echo "2" > modes.in')
    for m in mode.name.split("-"):
        m = split_digit_nondigit(m)
        mtype.extend(m[1].upper())
        modes.extend(m)
        os.system('echo "%03d %s %03d" >> modes.in' % (int(m[0]),
                                                       m[1].lower(),
                                                       int(m[2])))
    cdeg = min(max_cc_degrees(modes))

    # kernels generation
    for k in ["s", "p", "r"]:
        os.system('echo %s > input' % k)
        os.system('echo "%s" >> input' % cdeg)

        res = subprocess.Popen('%s < input' % cc_cstkernels, shell=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        output, error = res.communicate()

        if k == "s":
            kname = "vs"
        if k == "p":
            kname = "vp"
        if k == "r":
            kname = "rho"
        os.system('mv kernel.dat %s-kernel.dat' % kname)
        os.system('rm ast* cst.dat mcst.dat mdcpl.out raw.dat')

    # try CC with degree zero. e.g. code only works with 1s2-0s0 not 0s0-1s2
    if not os.path.exists("rho-kernel.dat"):
        for k in ["s", "p", "r"]:
            os.system('echo %s > input' % k)
            os.system('echo "%s" >> input' % 0)

            res = subprocess.Popen('%s < input' % cc_cstkernels, shell=True,
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.PIPE)
            output, error = res.communicate()

            if k == "s":
                kname = "vs"
            if k == "p":
                kname = "vp"
            if k == "r":
                kname = "rho"
            os.system('mv kernel.dat %s-kernel.dat' % kname)
            os.system('rm ast* cst.dat mcst.dat mdcpl.out raw.dat')

    os.system('rm input modes.in')
    if res.returncode != 0:
        raise Exception(error)
        print('error: %s') % res.returncode

    vp = []
    vs = []
    rho = []
    depth = []

    with open('rho-kernel.dat') as f:
        for line in f:
            cols = line.split()
            depth.append(float(cols[1]))
            rho.append(float(cols[0]))
    os.system('rm rho-kernel.dat')
    rhomax = np.amax(np.abs(rho))

    with open('vs-kernel.dat') as f:
        for line in f:
            cols = line.split()
            vs.append(float(cols[0]))
    os.system('rm vs-kernel.dat')
    vsmax = np.amax(np.abs(vs))
    intmax = np.amax([rhomax, vsmax])

    #if mode.type == 'S' and kernel in ('all', 'vp'):
    if mtype[0] != "T" and mtype[1] != "T":
        with open('vp-kernel.dat') as f:
            for line in f:
                cols = line.split()
                vp.append(float(cols[0]))
        vpmax = np.amax(np.abs(vp))
        intmax = np.amax([intmax, vpmax])
    os.system('rm vp-kernel.dat')

    ax.set_ylim(0, 6400)
    ax.set_xlim(-1.1*intmax, 1.1*intmax)

    ax.axhline(5700, color='k', lw=0.1)  # TZ
    ax.axhline(3480, color='k', lw=0.1)  # CMB
    ax.axhline(1220, color='k', lw=0.1)  # ICB
    ax.axvline(0, color='k', lw=0.1)

    if color == 'auto':
        rho_clr = 'grey'
        vs_clr = 'red'
        vp_clr = 'k'
    else:
        rho_clr = color[0]
        vs_clr = color[1]
        vp_clr = color[2]

    legend_ax = []
    legend = []
    if kernel in ('all', 'rho'):
        han1, = ax.plot(rho, depth, color=rho_clr, linewidth=linewidth)
        legend_ax.append(han1)
        legend.append(r'$\rho$')

    if kernel in ('all', 'vs'):
        han2, = ax.plot(vs, depth, color=vs_clr, linewidth=linewidth)
        legend_ax.append(han2)
        legend.append(r'$v_s$')

    #if mode.type == 'S' and kernel in ('all', 'vp'):
    if mtype[0] != "T" and mtype[1] != "T":
        if kernel in ('all', 'vp'):
            han3, = ax.plot(vp, depth, color=vp_clr, linewidth=linewidth)
            legend_ax.append(han3)
            legend.append(r'$v_p$')

    if 'ticks' in kwargs:
        ticks = kwargs['ticks']
    else:
        ticks = True

    if legend_show and ticks:
        if 'bbox_to_anchor' in kwargs:
            _bbox = kwargs['bbox_to_anchor']
        else:
            _bbox = (-0.08, -0.03)
        ax.legend(legend_ax, legend, frameon=False,
                  bbox_to_anchor=_bbox,
                  handlelength=0.5, handletextpad=0.1, loc='lower left',
                  fontsize=fontsize)

    elif legend_show and not ticks:
        if 'bbox_to_anchor' in kwargs:
            _bbox = kwargs['bbox_to_anchor']
        else:
            _bbox = (0.6, 1.2)
        ax.legend(legend_ax, legend, ncol=3,
                  handlelength=0.3, handletextpad=0.1, columnspacing=0.25,
                  loc='upper center', bbox_to_anchor=_bbox, frameon=False,
                  fontsize=fontsize, borderpad=0)
    if not ticks:
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        #ax.text(intmax*1.1, 5700, "TZ",  fontsize=fontsize*0.8, va="bottom", ha="right")
        #ax.text(intmax*1.1, 3480, "CMB", fontsize=fontsize*0.8, va="bottom", ha="right")
        #ax.text(intmax*1.1, 1220, "ICB", fontsize=fontsize*0.8, va="bottom", ha="right")
        ax.text(-intmax*1.1, 5700, "TZ",  fontsize=fontsize*0.8, va="bottom")
        ax.text(-intmax*1.1, 3480, "CMB", fontsize=fontsize*0.8, va="bottom")
        ax.text(-intmax*1.1, 1220, "ICB", fontsize=fontsize*0.8, va="bottom")
    if title:
        title = r"$_{%s}%s_{%s}$-$_{%s}%s_{%s}$" % (modes[0],modes[1],modes[2],
                                                    modes[3],modes[4],modes[5])
        ax.set_title(title, fontsize=fontsize*1.5, va="center")
    if show:
        plt.show()
    if savefig and ax is None:
        fig.savefig('%s_kernel.png' % mode.name, orientation='landscape',
                    dpi=400)

    os.chdir(cwd)
    # os.removedirs(tmp_path)
    os.system('rm -rf %s' % tmp_path)

    return fig, ax


def _get_plot_labels(label, colormap, colormap_index):
    marker = '^'
    linestyle = '-'
    if label == 'RR':
        marker = 'o'
        color = 'grey'
        linestyle = ':'

    elif label == 'REM':
        marker = 'o'
        color = 'grey'
        linestyle = '--'

    elif label == 'TZ':
        marker = 'x'
        color = 'grey'
        linestyle = '-.'

    elif label.startswith("S20RTS"):
        marker = '*'
        color = 'grey'
        linestyle = '-.'

    elif label.startswith("S40RTS"):
        marker = 's'
        color = 'grey'
        linestyle = ':'

    elif label.startswith("QRFSI12"):
        marker = 's'
        color = 'grey'
        linestyle = ':'

    elif label == 'AD':
        marker = 'd'
        color = 'grey'
        linestyle = '-'

    elif label == 'PK':
        marker = 'P'
        color = 'grey'
        linestyle = '--'
    else:
        if type(colormap) is str:
            color = colormap
        else:
            color = colormap[colormap_index]
            # color = next(colormap)

    return marker, color, linestyle


def _plot_map(clm, mode, kind, suptitle, html=False,
              kind_in_title=True, show_title=True, legend_show=True,
              ticks=False, weight="bold", llsvp=False, plates=False,
              outline_color="k", outline_thick=1, coastline_thick=1,
              meridians_thick=0.5, lon_0=180.0,
              **kwargs):

    bins = Bin()
    fig_config = {'figsize': (3, 1.8)}
    map_config = {'projection': 'kav7', 'lon_0': lon_0,
                  'resolution': 'c'}

    if 'smax' in kwargs:
        smax = kwargs['smax']

    if 'smin' in kwargs:
        smin = kwargs['smin']

    if 'fig' in kwargs:
        fig = kwargs['fig']
        fig.set_size_inches(8, 4)
    else:
        if html is False:
            fig = Figure(**fig_config)
        else:
            # print("html is ", html)
            fig = mpl.figure.Figure(**fig_config)

    if 'ax' in kwargs:
        ax = kwargs['ax']
    else:
        ax = None

    if 'show_colorbar' not in kwargs:
        show_colorbar = True
    else:
        show_colorbar = False

    if 'cmap' in kwargs:
        cmap = kwargs['cmap']
        if type(cmap) == str:
            if hasattr(cm, cmap):
                cmap = getattr(cm, cmap)
            elif cmap.lower() == 'arwen' or cmap.lower() == 'yr':
                colors = [(1, 0, 0), (1, 1, 0), (0, 0, 1)]
                cmap_name = 'sf_cmap'
                cmap = LinearSegmentedColormap.from_list(cmap_name,
                                                         colors, N=12)
            else:
                colors = [(0.8, 0, 0.2), (1, 0, 0), (1, 0.6, 0.1), # DR, R, LOr
                          (1, 1, 0),# Y
                          (0.3, 0.6, 1), (0.3, 0.3, 0.9), (0.2, 0, 0.7)] # LB, MB, DB
                #colors = [(1, 0, 0), (1, 1, 0), (0, 0, 1)] # R -> Y -> B
                if kind == "dst":
                    colors = colors[::-1]  # B -> Y -> R
                cmap_name = 'sf_cmap'
                cmap = LinearSegmentedColormap.from_list(cmap_name,
                                                         colors, N=12)
        else:
            # colors = [(1, 0, 0), (0.878, 1, 1), (0, 0, 1)]
            colors = cmap
            cmap_name = 'sf_cmap'
            cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=12)
    else:
        colors = [(1, 0, 0), (0.9, 1, 0.9), (0, 0, 1)]  # R -> Y -> B
        if kind == "dst":
            colors = colors[::-1]  # B -> Y -> R

        cmap_name = 'sf_cmap'
        cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=12)

    for arg, val in kwargs.items():
        if arg == 'kind_in_title':
            kind_in_title = False

        if arg == 'show_title':
            show_title = False

        if arg == 'weight':  # letter weight
            weight = val

        if arg == 'llsvp':  # llsvp outline
            llsvp = val

        if arg == 'plates':  # NEW plate tectonics outline
            plates = val

        if arg == 'outline_color':  # llsvp/plates outline color
            outline_color = val

        if arg == 'outline_thick':  # llsvp/plates outline thickness
            outline_thick = val

        if arg == 'coastline_thick':  # coastline thickness
            coastline_thick = val

        if arg == 'meridians_thick': # meridian thickness
            meridians_thick = val

    lons = np.arange(-180, 180)
    lats = np.arange(-90, 90)
    g = np.zeros([lats.size, lons.size])

    for j, lon in enumerate(lons):
        for i, lat in enumerate(lats):
            if lon < 0:
                lon = 360. + lon
            x = clm.expand(lat=[lat], lon=[lon])
            g[i, j] = x[0]

    vminp = np.floor(g.min()*10.)/10.
    vmaxp = np.floor(g.max()*10.)/10.

    if 'vmin' in kwargs:
        if kwargs['vmin'] is not None:
            vminp = kwargs['vmin']
    if 'vmax' in kwargs:
        if kwargs['vmax'] is not None:
            vmaxp = kwargs['vmax']

    if vminp < 0 and vmaxp > 0:
        mpnt = 0.
    else:
        mpnt = 0.5*(vminp + vmaxp)
        # mpnt = 1 - vmaxp/(vmaxp + abs(vminp))

    cmap = shiftedColorMap(cmap, start=0,
                           midpoint=(mpnt-vminp)/(vmaxp-vminp),
                           stop=1.0, name='shiftedcmap')

    if '-' not in mode.name:
        if ax is None:
            if os.path.exists(bins.sc_cstkernels):
                gs = gridspec.GridSpec(1, 2, width_ratios=[0.8, 3.2],
                                       wspace=0.05)

                ax0 = fig.add_subplot(gs[0])
                ax = fig.add_subplot(gs[1])
                sens_kernel(mode, title=False, ax=ax0, kind=kind, **kwargs)

                if show_colorbar is True:
                    ax_cb = fig.add_axes([0.4, 0.1, 0.4, 0.04])

            else:
                ax = fig.add_subplot()

                if show_colorbar is True:
                    # To do: Have to change the values here to make it centered
                    ax_cb = fig.add_axes([0.35, 0.1, 0.3, 0.02])

        elif show_colorbar is True:
            # ax_cb = fig.add_axes([0.3, 0.1, 0.3, 0.02])
            ax_cb = fig.add_axes([0.25, 0.1, 0.4, 0.04])

        m = Basemap(ax=ax, **map_config)
        if show_colorbar is True:
            if smin == smax:
                _sdegs = "$s$=${}$".format(smin)
            else:
                _sdegs = "$s$=${}$-${}$".format(smin, smax)

            if kind_in_title:
                title = "%s %s \n $_{%s}%s_{%s},$ %s"
                title = title % (kind, suptitle, mode.n, mode.type,
                                 mode.l, _sdegs)
            elif not show_title:
                title = "$_{%s}%s_{%s},$ %s"
                title = title % (mode.n, mode.type, mode.l, _sdegs)
            else:
                title = "%s \n $_{%s}%s_{%s},$ %s"
                title = title % (suptitle, mode.n, mode.type,
                                 mode.l,_sdegs)
        else:
            if kind_in_title:
                title = "%s %s \n $_{%s}%s_{%s},$ $s_{max}=%s$"
                title = title % (kind, suptitle, mode.n, mode.type,
                                 mode.l, smax)
            elif not show_title:
                title = "$_{%s}%s_{%s},$ $s_{max}=%s$"
                title = title % (mode.n, mode.type, mode.l, smax)
            else:
                title = "%s \n $_{%s}%s_{%s},$ $s_{max}=%s$"
                title = title % (suptitle, mode.n, mode.type, mode.l, smax)
        if 'filename' in kwargs:
            ax.set_title('%s \n %s' % (kwargs['filename'], title))
        elif not kind_in_title:
            ax.set_title(title, weight="bold")
        else:
            ax.set_title(title)

    else: # CC splitting func
        if ax is None:
            if os.path.exists(bins.sc_cstkernels):
                gs = gridspec.GridSpec(1, 2, width_ratios=[0.8, 3.2],
                                       wspace=0.05)

                ax0 = fig.add_subplot(gs[0])
                ax = fig.add_subplot(gs[1])
                sensC_kernel(mode, title=False, ax=ax0, **kwargs)

                if show_colorbar is True:
                    ax_cb = fig.add_axes([0.4, 0.1, 0.4, 0.04])

            else:
                ax = fig.add_subplot()

                if show_colorbar is True:
                    # To do: Have to change the values here to make it centered
                    ax_cb = fig.add_axes([0.35, 0.1, 0.3, 0.02])

        m = Basemap(ax=ax, **map_config)
        ccn = split_digit_nondigit(mode.name)
        if smin == smax:
            _sdegs = "$s$=${}$".format(smin)
        else:
            _sdegs = "$s$=${}$-${}$".format(smin, smax)
        if kind_in_title:
            title = "%s %s\n$_{%s}%s_{%s}$-$_{%s}%s_{%s},$ %s"
            title = title % (kind, suptitle,
                             int(ccn[0]), str(ccn[1]).upper(), int(ccn[2]),
                             int(ccn[4]), str(ccn[5]).upper(), int(ccn[6]),
                             _sdegs)
        elif not show_title:
            title = "$_{%s}%s_{%s}$-$_{%s}%s_{%s},$ %s"
            title = title % (int(ccn[0]), str(ccn[1]).upper(), int(ccn[2]),
                             int(ccn[4]), str(ccn[5]).upper(), int(ccn[6]),
                             _sdegs)
        else:
            title = "%s\n$_{%s}%s_{%s}$-$_{%s}%s_{%s},$ %s"
            title = title % (suptitle,
                             int(ccn[0]), str(ccn[1]).upper(), int(ccn[2]),
                             int(ccn[4]), str(ccn[5]).upper(), int(ccn[6]),
                             _sdegs)
        if 'filename' in kwargs:
            ax.set_title('%s \n %s' % (kwargs['filename'], title))
        elif not kind_in_title:
            ax.set_title(title, weight="bold")
        else:
            ax.set_title(title)

    cp = m.drawmapboundary()
    m.drawcoastlines(linewidth=coastline_thick)
    m.drawparallels(np.arange(-90., 120., 60.), linewidth=meridians_thick, zorder=1, dashes=(None,None))
    m.drawmeridians(np.arange(0., 420., 60.), linewidth=meridians_thick, zorder=1, dashes=(None,None))

    if plates: # scattering has no wrapping depending on lon_0
        path = tplates.__path__[0]
        f = '%s/original/PB2002_plates.dig.ALL.txt' % path
        f = open(f, 'r')
        tp = [[l for l in line.split(",")] for line in f.readlines() if not line.startswith("#")]
        tp = [list(i) for i in zip(*tp)]
        tp = [[float(line) for line in token] for token in tp[0:2]]
        x, y = m(tp[0], tp[1])
        ax.scatter(x, y, facecolors=outline_color, linewidths=0, s=outline_thick,
                   marker=",", zorder=500)
        #ax.scatter(x, y, c=outline_color, s=outline_thick, zorder=500)

    if llsvp: # scattering has no wrapping depending on lon_0
        path = LLSVP.__path__[0]
        f = '%s/LLSVP_countour.dat' % path
        f = open(f, 'r')
        tp = [[float(l) for l in line.split()] for line in f.readlines()]
        tp = [list(i) for i in zip(*tp)]
        x, y = m(tp[0], tp[1])
        #ax.scatter(x, y, c=outline_color, s=outline_thick, zorder=500)
        ax.scatter(x, y, facecolors=outline_color, linewidths=0, s=outline_thick,
                   marker=",", zorder=500)
    if 'vmin' in kwargs or 'vmax' in kwargs:
        s = np.clip(g, vminp, vmaxp)
        s = m.transform_scalar(s, lons, lats, 1000, 500)
    else:
        s = m.transform_scalar(g, lons, lats, 1000, 500)

    im = m.imshow(s, vmin=vminp, vmax=vmaxp, cmap=cmap, clip_path=cp)

    if vminp < 0 and vmaxp > 0:
        ticks = [vminp, 0., vmaxp]
    else:
        ticks = [vminp, vmaxp]

    # plot it and set up / format labels
    if show_colorbar is True:
        cb = fig.colorbar(im, cax=ax_cb, ticks=ticks, format='%3.1f',
                          orientation='horizontal', extend='both')
        if html is True:
            pass
            # cb = plt.colorbar(im, cax=ax_cb, ticks=ticks, format='%3.1f',
            #                   orientation='horizontal', extend='both')
            # cb.clim(s.min(), s.max())
        else:

            cb.set_clim(s.min(), s.max())

        cb.ax.set_title(r'$\mu$Hz', x=1.2, y=-0.7)

    if 'savefig' in kwargs:
        save = kwargs['savefig']
        if 'filename' in kwargs and save:
            name = kwargs['filename']
            fname = '%s_%s_%s' % (mode.name, kind, name)
            fig.savefig('%s.png' % fname, orientation='landscape', dpi=400,
                        bbox_inches="tight", pad_inches=0.01, transparent=True)
        elif save:
            fname = '%s_%s' % (mode.name, kind)
            fig.savefig('%s.png' % fname, orientation='landscape', dpi=400,
                        bbox_inches="tight", pad_inches=0.01, transparent=True)

    return im, fig


def _plot_coeffs(coeffs, errors, mode_name, label, modes, kind,
                 smin=-1, smax=-1, colormap=None, colormap_index=None,
                 d00=None, d00_err=None, html=False, **kwargs):

    """
    Setup Figures, Initiate Axes
    accepted values for kwargs:
    label
    fig_size: list of 2 integers
    error_bars
    """
    if 'error_bars' in kwargs:
        error_bars = kwargs['error_bars']
    else:
        error_bars = 'symmetric'
    # print('html is', html)
    a_i = 0
    Ni = 0
    Nf = len(coeffs)
    sdegs = sorted([int(x) for x in coeffs.keys()])
    if smin != -1 and smax != -1:
        if smin in sdegs:
            Ni = sdegs.index(smin)
        if smax <= max(sdegs):
            Nf = sdegs.index(smax) + 1

    elif smax != -1 and smax <= max(sdegs):
        Nf = sdegs.index(smax) + 1

    elif smin in sdegs:
        Ni = sdegs.index(smin)

    N = Nf - Ni

    if mode_name not in modes:
        modes[mode_name] = AttribDict()
        # if html is False:
        #     fig = plt.figure()
        # else:
        fig = mpl.figure.Figure()
        gs = gridspec.GridSpec(N, 1)

        # Try to get mode name
        try:
            mN, mK, mL = split_digit_nondigit(mode_name)
            _title = r"%s $_{%g}%s_{%g}$" % (kind, int(mN), mK, int(mL))
        except Exception:
            _title = r"%s $%s$" % (kind, mode_name)

        if 'title' in kwargs:
            _title += "\n%s" % kwargs['title']

        fig.suptitle(_title)
        modes[mode_name]['fig'] = fig
        modes[mode_name]['axes'] = AttribDict()

    else:
        fig = modes[mode_name]['fig']

    axes = modes[mode_name]['axes']

    # init plots
    marker, c, ls = _get_plot_labels(label, colormap, colormap_index)

    """
    Setup Axes
    """
    # Create axes dictionary here, with ax entry for each degree s.
    # degree is the current s value, cval are the corresponding splitting
    # function coefficients (cst or dst)
    coeffslist = sort_human(list(coeffs.keys()))
    for degree in coeffslist[Ni:Nf]:
        cval = coeffs[degree]
        # Check if there exists a ax for the current s, else create oone
        if degree not in axes:
            axes[degree] = AttribDict()
            axes[degree]['xmax'] = 0
            axes[degree]['ymax'] = 0
            # if html is False:
            #     axes[degree]['ax'] = plt.subplot2grid((N, 2), (a_i, 0),
            #                                           colspan=2)
            # else:
            axes[degree]['ax'] = fig.add_subplot(gs[a_i, :])
            axes[degree]['ax'].set_title(r"$s=%i$" % int(degree))
            a_i += 1

        # Get plot limits here
        if abs(cval).max() > abs(axes[degree]['ymax']):
            axes[degree]['ymax'] = abs(cval).max()

        if len(cval) > axes[degree]['xmax']:
            xlim = len(cval) - 0.5
            if degree == '0':
                xlim += 1
            axes[degree]['xmax'] = xlim

    """
    Plotting loop
    """
    for degree in coeffslist[Ni:Nf]:
        cval = coeffs[degree]
        try:
            if error_bars == 'symmetric':
                eval = errors[degree]['uncertainty']
            else:
                e1 = np.absolute(errors[degree]['lower_uncertainty'])
                e2 = np.absolute(errors[degree]['upper_uncertainty'])
                eval = np.array([e1, e2])
        except Exception:
            eval = [None] * len(cval)

        if degree == '0':
            if d00 is not None:
                cval = np.append(cval, d00)
                try:
                    if error_bars == 'symmetric':
                        eval = np.append(eval, d00_err['uncertainty'])
                    else:
                        e1 = np.absolute(d00_err['lower_uncertainty'])
                        e2 = np.absolute(d00_err['upper_uncertainty'])
                        eval_n = np.array([e1, e2])
                        eval = np.append(eval, [eval_n[0], eval_n[1]], axis=1)
                except Exception:
                    eval = np.append(eval, None)

        x = np.arange(len(cval))
        ax = axes[degree].ax
        xlim = axes[degree].xmax

        ax.hlines(0, -0.5, xlim, linestyle='-', color='black',
                  linewidth=.45)

        if label is not None:
            kwargs['label'] = label

        if degree == '0':
            ax.set_xlim((-0.5, xlim))
            ax.hlines(0, -0.5, xlim, linestyle='-', color='black',
                      linewidth=.45)
        else:
            ax.set_xlim((-0.5, xlim))
            ax.hlines(0, -0.5, xlim, linestyle='-', color='black',
                      linewidth=.45)

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.errorbar(x, cval, yerr=eval, capsize=3, marker=marker,
                    linestyle=ls, color=c, linewidth=0.8, label=label)
    if 'label' in kwargs:
        if 'bbox_to_anchor' in kwargs:
            bbox = kwargs['bbox_to_anchor']
        else:
            bbox = (1.05, 1)
        if 'loc' in kwargs:
            loc = kwargs['loc']
        else:
            loc = 2
        # removing whiskers from legend
        ehandles, elabels = ax.get_legend_handles_labels()
        ehandles = [h[0] for h in ehandles]
        degs = [int(y) for y in coeffslist[Ni:Nf]]
        ax = axes[str(min(degs))].ax
        ax.legend(ehandles, elabels,
                  bbox_to_anchor=bbox, loc=loc, borderaxespad=0.)

    if 'fig_size' in kwargs:
        fig.set_size_inches(kwargs['fig_size'])

    return modes
