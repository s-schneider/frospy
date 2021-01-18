import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import matplotlib.markers as mmarkers

import os
import re
from os.path import join
from obspy.core import AttribDict
from collections import OrderedDict
import time
import pandas as pd

from frospy.postprocessing.read import (query_ev_seg_misfits,
                                      query_inversion_summary,
                                      read_cst_files)
from frospy.core.splittingfunc import loadmodel, Set, SplittingFunc
from frospy.core.splittingfunc.plot import sens_kernel
from frospy.core.database.query import get_db_tables, db_query
from frospy.core.modes import format_name
from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes, Mode
from frospy.core.segment import Segment
from frospy.core.segment import read as read_seg
from frospy.core.setup.settings import read as read_setup
from frospy.core.setup.settings import Setup

from frospy.util.read import read_std_cat, read_std_inv
from frospy.util.base import (format_list, update_progress, split_digit_nondigit,
                            is_number, is_mode, uniq_modes, sort_catalog,
                            max_cc_degrees, cst2fQ)
from frospy.plot.nmplt import get_iter_colormap
import glob
import obspy


def get_db_path(run_dir):
    if type(run_dir) == list:
        if not os.path.exists(run_dir[0]):
            print('run_dir does not exist!')
            return
        if run_dir[1] == 1:
            inv_folder = "inversion_out"
        else:
            inv_folder = "inversion_out_%s" % run_dir[1]

        # Overwriting run with the folder_name
        run_dir = run_dir[0]
        db_path = join(run_dir, inv_folder, 'inversion.db')
        if not os.path.exists(db_path):
            db_path = None

    else:
        if not os.path.exists(run_dir):
            print('run_dir does not exist!')
            return None
        db_paths = glob.glob(join(run_dir, 'inversion_out*'))

        if len(db_paths) == 0:
            print('inversion_out directory does not exist!')
            return None

        db_path = join(run_dir, 'inversion_out', 'inversion.db')

        if not os.path.exists(db_path):
            print('inversion.db file not found')
            return None

    return db_path


def inv_info(run_dir, fig_size=(12, 7), font_size=9, table_scale='auto',
             plot_kernel=True, verbose=False):
    # Read the setup.pickle file here and get the following info:
    #  - modes
    #  - degrees

    runs = glob.glob('{path}/*'.format(path=run_dir))
    runs.sort()

    if len(runs) == 0:
        print('No runs found in %s' % run_dir)
        return

    runs_old = glob.glob('{path}/old/*'.format(path=run_dir))
    if verbose:
        msg = '%i runs in old' % len(runs_old)
        print(msg)

    cell_text = []
    columns = ('Run', '$N_{events}$', '$N_{T}$', '$N_{R}$', '$N_{Z}$',
               '$N_{params}$', 'Modes, $l_{max}$', 'Startmodel', '#')
    colWidths = [0.12, 0.05, 0.05, 0.05, 0.05, 0.1, 0.15, 0.1, 0.05]

    modes_kernel = Modes()
    for run in runs:
        if verbose:
            print('Run:\t\t%s' % run)
        try:
            setup = read_setup('{path}/setup.pickle'.format(path=run))
        except Exception as e:
            if verbose is True:
                print('  Error\t%s' % e)

        if setup is None:
            continue
        modes_kernel += setup.modes
        try:
            run_name = os.path.basename(run)
            r_final = glob.glob('{path}/inversion_out'.format(path=run))
            r = glob.glob('{path}/inversion_out_*'.format(path=run))
            r = r_final + r
            r = r[::-1]
            for i, inv in enumerate(r):
                if inv.endswith('out'):
                    index = 1
                else:
                    index = int(inv.split('inversion_out_')[-1])
                if verbose:
                    print('  Checking\t%s' % inv.split('/')[-1])
                db_path = join(inv, 'inversion.db')
                if len(setup.modes) == 1:
                    modes_d = ''
                else:
                    modes_d = '\n'

                NRec = {'Z': 0, 'R': 0, 'T': 0}
                if not os.path.exists(db_path):
                    Nev = '-'
                else:
                    tables = get_db_tables(db_path)
                    ref_damp = float([x for x in tables if is_number(x)][0])
                    for t in tables:
                        if t == 'initial' or is_number(t):
                            continue

                        sel = {'iter': '1', 'damp': ref_damp}
                        events = db_query(db_path, t, select='event',
                                          condition_dict=sel)
                        Nev = len(set(events))

                        if is_mode(t) is not False or t == 'Set':
                            # Cross-Coupling behavior adding
                            if is_mode(t) == 'S' or setup.modes[0].type == 'S':
                                NRec['Z'] = len(events)
                            else:
                                NRec['T'] = len(events)

                        elif t[0] in ('T', 'Z', 'R'):
                            NRec[t[0]] += len(events)
                    if len(setup.modes_sc) > 1:
                        for m, d in setup.modes_sc.iteritems():
                            modes_d += '%s, %i \n ' % (m, d)
                    else:
                        modes_d += '%s, %i' % (setup.modes_sc.items()[0])

                    if setup.modes_cc is not None:
                        for m, d in setup.modes_cc.iteritems():
                            modes_d += '%s, %i-%i \n' % (m, d[0], d[1])

                if setup.startmodel in ('S20RTS', 'PREM'):
                    startmodel = setup.startmodel
                else:
                    startmodel = 'custom'

                cell_text.append([run_name, Nev,
                                  NRec['T'], NRec['R'],  NRec['Z'],
                                  setup.model_size, modes_d, startmodel,
                                  index])
        except Exception as e:
            if verbose:
                print("Error \n%s\n in %s" % (e, run))
            else:
                pass

    modes_kernel = uniq_modes(modes_kernel)

    df = pd.DataFrame(cell_text, columns=columns)

    fig = plt.figure()
    fig.set_size_inches(fig_size)

    gs = gridspec.GridSpec(1, len(modes_kernel)+5)
    if plot_kernel is True:
        ax_kernel = []
        for i, m in enumerate(modes_kernel):
            ax_kernel.append(plt.subplot(gs[0, i]))
        for i, m in enumerate(modes_kernel):
            ticks = legend_show = True
            if i > 0:
                ticks = legend_show = False
            sens_kernel(m, ax=ax_kernel[i], ticks=ticks,
                        legend_show=legend_show)

    ax_table = plt.subplot(gs[0, len(modes_kernel):])
    ax_table.axis('off')
    ax_table.axis('tight')
    table = ax_table.table(cellText=df.values, colLabels=df.columns,
                           loc='center', colWidths=colWidths)

    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    if table_scale == 'auto':
        if len(modes_kernel) == 1:
            table.scale(1.5, 1.5)
        else:
            table.scale(1.5, len(modes_kernel) * 2.1)
    else:
        table.scale(table_scale[0], table_scale[1])

    for run in runs:
        if verbose:
            print('\n##### Generating Summary #####\n\n')
        job_files = glob.glob('%s/main*err' % run)
        job_ids = " |"
        for j in job_files:
            job_ids += " %s |" % j.split('main')[-1].split('.')[1]
        print("%s %s" % ('/'.join(run.split('/')[-3:]), job_ids))

        try:
            inv_summary(None, run_dir=run, damping='all',
                        print_summary=True, plot=False,
                        verbose=verbose)
        except Exception as e:
            print('Error: %s' % e)

    return


def inv_summary(schema, setup=None, run_dir=None, iteration='all',
                damping='all', plot=True, title=None, histtype='bar',
                bins=None, verbose=False, lcurve=True, fig_size=(10, 7),
                print_summary=False, plot_events=False, plot_stations=False,
                **plotargs):

    if setup is None and run_dir is None:
        raise IOError('No setup or run_dir provided')
    if setup is not None:
        run_dir = join(setup.rundir, setup.inversion_outdir)

    db_path = get_db_path(run_dir)
    if db_path is None:
        return

    if type(schema) is str:
        schema = [schema]
    if verbose:
        t1 = time.time()
    # read inversion.db
    avmisfit, Msize, effEV, eps, N, dampings = query_inversion_summary(db_path,
                                                                       verbose)
    if verbose:
        print("\nNumber of files per Damping:")
        for _d, _val in N.iteritems():
            print("\t%s: %s" % (_d, _val))

    damp_keys = N.keys()
    damp_keys.remove('initial')
    iterations = N[damp_keys[0]].keys()
    if verbose:
        print("\nIterations: %s \n" % iterations)
    N_r = N['initial']  # [max(iterations)]

    # printing misfits
    # del avmisfit['initial']
    if plot is not True or print_summary is True:
        msg = "| damp\t\t| init\t\t| final\t\t| N \t| Iterations"
        print(msg)
        # for damp, iters in avmisfit.iteritems():
        for damp in dampings[:-1]:
            # print(damp)
            if damping != 'all':
                if type(damping) is not list:
                    if float(damp) != float(damping):
                        continue
                else:
                    if float(damp) not in [float(x) for x in damping]:
                        continue
            enditer = max(avmisfit[damp].keys())

            msg = "|%8.4f\t|%8.4f\t|%8.4f\t| %i\t| %i\t"
            print(msg % (
                damp, avmisfit['initial'], avmisfit[damp][enditer],
                N_r, enditer))
        if plot is not True:
            return

    if verbose:
        t2 = time.time()
        print("Time for database search: %f s" % abs(t1-t2))

    fig, ax = plt.subplots()

    # Number of rows in plot
    N_rows = 3 * len(schema) + 3

    ax = AttribDict()
    ax['Ave_mf'] = plt.subplot2grid((N_rows, 8), (0, 0), rowspan=2, colspan=2)
    ax['Msize'] = plt.subplot2grid((N_rows, 8), (0, 3), rowspan=2, colspan=2)
    ax['EffEV'] = plt.subplot2grid((N_rows, 8), (0, 6), rowspan=2, colspan=2)
    ax['Src_hist'] = {}
    ax['Src_mf'] = {}
    ax['mf_hist'] = {}

    if verbose:
        t3 = time.time()
        print("Time for Figure creation: %f s" % abs(t3-t2))

    xlabel = 'iteration'
    ylabel = r'$\frac{\sum mf^2}{\sum data^2}$'
    plot_graph(ax.Ave_mf, avmisfit, damping, iteration, 'Average Misfit',
               xlabel, ylabel, legend=False)

    if damping != 'all' or lcurve is not True:
        ylabel = r'$\sum cst_i^2$'
        plot_graph(ax.Msize, Msize, damping, iteration, 'Modelsize',
                   xlabel, ylabel, legend=False)
    else:
        plot_L_curve(ax.Msize, avmisfit, Msize, iterations, iteration)

    ylabel = r'$\sum diag(ResMat)$'
    plot_graph(ax.EffEV, effEV, damping, iteration, 'Eff. EV',
               xlabel, ylabel)

    # get number of files
    N_records = {}
    for i, s in enumerate(schema):
        # Query Inv(ersion) object
        if verbose:
            print("Searching for: \n%s\n%s\n%s" % (s, iterations, dampings))
        Inv = query_ev_seg_misfits(db_path, s, iterations, dampings)

        N_records[s] = 0
        for val in Inv.stats.itervalues():
            N_records[s] += val['n']
        ax['Src_hist'][s] = plt.subplot2grid((N_rows, 8), (i*3 + 3, 0),
                                             rowspan=2, colspan=2)
        ax['Src_mf'][s] = plt.subplot2grid((N_rows, 8), (i*3 + 3, 3),
                                           rowspan=2, colspan=2)
        ax['mf_hist'][s] = plt.subplot2grid((N_rows, 8), (i*3 + 3, 6),
                                            rowspan=2, colspan=2)

        if i == 0:
            legend = True
        else:
            legend = False

        plot_N_per_event(ax.Src_hist[s], Inv,
                         ylabel='No. of Records %s' % s, N=i)

        plot_mf_per_event(ax.Src_mf[s], Inv, damping=damping,
                          N=i, **plotargs)

        plot_histogram(ax.mf_hist[s], Inv, 'Misfit', 'No. of Records',
                       bins=bins, histtype=histtype, legend=legend,
                       damping=damping, N=i, title='mean', **plotargs)

    if title:
        title = '%s, using Records - ' % (title)
        for s, n in N_records.iteritems():
            title += ' %s: %i/%i  ' % (s, n, N_r)
    else:
        title = 'Records - '
        for s, n in N_records.iteritems():
            title += ' %s: %i/%i  ' % (s, n, N_r)

    fig.suptitle(title)
    fig.set_size_inches(fig_size)

    if verbose:
        t4 = time.time()
        print("Time for Plotting: %f s" % abs(t3-t4))

    if plot_stations is True:
        setup = read_setup(join(run_dir, 'setup.pickle'))
        inv = read_std_inv()
        inv_plot = obspy.Inventory(networks=[], source="nmpy")

        seg = Segment()
        for dir in setup.segmentsdir.iterkeys():
            for sdir in glob.glob("%s/???????" % join(run_dir, dir)):
                for segments in glob.glob("%s/*%s" % (sdir,
                                                      setup.seg_suffix[dir])):
                    seg += read_seg(segments)
        for p in seg:
            inv_plot += inv.select(station=p.stats.station)
        inv_plot.plot()

    plt.show()

    if plot_events is True:
        if damping != 'all':
            dampings = format_list(damping, output='float')
        else:
            dampings.remove('initial')

        cat = read_std_cat(Inv.stats.keys())
        cat = sort_catalog(cat, 'time')
        cat.plot()
        msg = '  id    |      origin       | lat, lon | magnitude | depth '
        msg += '| N of events | Mf \n'
        for event in cat:
            cmt = event.event_descriptions[-1].text
            event_stats = Inv.stats[cmt]

            ev_str = [event.__str__().split()[1][0:19]]
            ev_str += event.__str__().split()[2:8]
            msg += '%s | ' % cmt
            msg += ' '.join(ev_str)
            msg += ' | %.1f km' % (event.origins[0].depth / 1000.)
            msg += ' | %i' % event_stats.n
            msg += ' | %.2f' % event_stats['misfit']['initial']
            for d in dampings:
                msg += ' | %.2f' % event_stats['misfit'][d]
            msg += '\n'
        print(msg)

    return fig, ax


def plot_L_curve(ax, mf, msize, iterations, iter_no):

    if iter_no != 'all':
        it = int(iter_no[-1])
        title = "L-Curve, Iteration: %s" % it
    else:
        title = "L-Curve"

    damps = mf.keys()
    damps.remove('initial')
    damps.sort()
    # damps = np.array(damps)

    L = []
    for d in damps:
        if iter_no == 'all':
            it = mf[d].keys()[-1]
        L.append([mf[d][it], msize[d][it]])
        x = mf[d][it]
        y = msize[d][it]
        ax.text(x, y, d,
                horizontalalignment='center',
                verticalalignment='bottom')

    L = np.array(L)
    ax.plot(L.transpose()[0], L.transpose()[1], color='blue')
    sc = ax.scatter(L.transpose()[0], L.transpose()[1], s=100, marker='v')
    ax.set_xlabel(r'$\frac{\sum mf^2}{\sum data^2}$')
    ax.set_ylabel(r'$\sum cst_i^2$')
    ax.grid(True)

    ax.set_title(title)
    fig = plt.gcf()

    names = [str(d) for d in damps]
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "d: {}".format(" ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(damps[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    return


def plot_N_per_event(ax, Inv, title='Record Distribution',
                     xlabel='Events', ylabel='No. of Records',
                     rotate_xlabel=90, N=0, **args):
    if 'alpha' not in args:
        alpha = 1
    else:
        alpha = 1

    xtick_label = []
    y = []

    # Calcualte Mean Number of events
    N_mean = 0
    for event in Inv.stats.itervalues():
        N_mean += event.n
    N_mean = N_mean / float(len(Inv.stats))

    for cmt, values in Inv.stats.iteritems():
        if values.n > 2 * N_mean or len(Inv.stats.keys()) <= 10:
            xtick_label += [str(cmt)]
        else:
            xtick_label += ['']

        y += [values['n']]
    x = range(len(y))

    ax.bar(x, y, alpha)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    if ylabel is not None:
        ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(y)))
    ax.set_xticklabels(xtick_label)

    for tick in ax.get_xticklabels():
        if rotate_xlabel is not False:
            tick.set_rotation(rotate_xlabel)
    if N == 0 and title is not None:
        ax.set_title(title)
    return


def plot_mf_per_event(ax, Inv, title='Events Misfit',
                      xlabel='Events', ylabel='misfit', yvalue='misfit',
                      rotate_xlabel=90, damping='all', N=0, **args):
    if 'alpha' not in args:
        args['alpha'] = 0.8

    x_ticklabel = []
    y = []
    x_init = []
    y_init = []

    keys = Inv.stats.values()[0]['misfit'].keys()
    keys.sort()
    keys.reverse()

    if damping != 'all':
        damping = format_list(damping, output='float')
    else:
        damping = keys[:]
        damping.remove('initial')
    step = len(damping)

    cmts = Inv.stats.keys()
    cmts.sort()

    # Calcualte Mean misfit of events
    mf_mean = 0
    for cmt, values in Inv.stats.iteritems():
        mf_mean += values['misfit']['initial']
    mf_mean = mf_mean / float(len(Inv.stats))

    x_n = 0
    for cmt, values in Inv.stats.iteritems():
        for k in keys:
            if k == 'initial':
                x_init += [x_n]
                y_init += [values['misfit'][k]]
                if (
                    values['misfit'][k] > 2 * mf_mean or
                    len(Inv.stats.keys()) <= 10
                   ):
                    x_ticklabel += [str(cmt)]
                else:
                    x_ticklabel += ['']
            elif k in damping:
                y += [values['misfit'][k]]
        x_n += 1

    ax.bar(x_init, y_init, color='black')
    for i, k in enumerate(damping):
        ax.bar(x_init, y[i::step], **args)

    ax.set_xticks(x_init)
    ax.set_xticklabels(x_ticklabel)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    if ylabel is not None:
        ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if rotate_xlabel is not False:
        for tick in ax.get_xticklabels():
            tick.set_rotation(rotate_xlabel)
    if N == 0 and title is not None:
        ax.set_title(title)

    return


def plot_bar(ax, event, title=None, xlabel=None, ylabel=None, yvalue='N',
             rotate_xlabel=False, **args):
    x = []
    y = []
    for cmt, values in event.iteritems():
        x += [str(cmt)]

        if yvalue == 'N':
            y += [len(values.values()[0])]

        elif yvalue == 'misfit':
            N = len(values[max(values.keys())])
            mf = sum(values[max(values.keys())])
            y += [mf / float(N)]

    ax.bar(x, y, **args)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    if ylabel is not None:
        ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if rotate_xlabel is not False:
        for tick in ax.get_xticklabels():
            tick.set_rotation(rotate_xlabel)
    if title is not None:
        ax.set_title(title)

    return


def plot_histogram(ax, Inv, xlabel=None, ylabel=None, legend=None,
                   damping='all', title=None, histtype='bar', bins=None, N=0,
                   **args):
    if 'alpha' not in args:
        args['alpha'] = 0.8

    if damping != 'all':
        damping = format_list(damping, 'float', output='float')
        damping += ['initial']

    N = len(Inv.allmisfits.initial)
    if bins is None:
        bins = int(N * 0.5)
    if bins < 20:
        bins = 20

    keys = Inv.allmisfits.keys()
    keys.sort()
    keys.reverse()
    sorted(keys, key=str)

    values = np.array([])
    for key in keys:
        if damping != 'all':
            if key not in damping:
                continue

        Inv.allmisfits[key]

        if histtype == 'line':
            y, binEdges = np.histogram(Inv.allmisfits[key], bins)
            bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
            ax.plot(bincenters, y, **args)
        elif histtype == 'bar':
            if key == 'initial':
                ax.hist(Inv.allmisfits[key], bins=bins, histtype=histtype,
                        label=str(key), color='black', **args)
            else:
                ax.hist(Inv.allmisfits[key], bins=bins, histtype=histtype,
                        label='d=%s' % str(key), **args)
        values = np.hstack((values, Inv.allmisfits[key]))
#        if key == 'initial':
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if N == 0 and title is not None:
        ax.set_title(title)
    elif title == 'mean':
        if len(damping) == 2:
            title = "Mean: %.2f" % np.mean(values)
            ax.set_title(title)

    ax.grid(True)

    return


def plot_graph(ax, Dict, damping, iteration, title, xlabel, ylabel,
               marker='v', legend=True):

    if damping == 'all':
        damping = Dict.keys()
        damping.sort()
        damping.reverse()
    else:
        damping = format_list(damping, 'float', output='float')
        damping += ['initial']

    xticks = []
    for d in damping:
        if d == 'initial':
            continue
        try:
            values = [Dict['initial']]
            iters = [0]
        except KeyError:
            values = []
            iters = []

        if iteration == 'all':
            iter_no = Dict[d].keys()
            for i in iter_no:
                if i not in xticks:
                    xticks.append(i)

        for i in iter_no:
            values += [Dict[d][i]]
            iters += [i]

        # label format
        if (d > 1e-5 and d < 10000) or d == 0:
            dformat = "d=%g"
        else:
            dformat = "d=%4.2e"
        ax.plot(iters, values, marker=marker, label=dformat % d)
    xticks.sort()
    ax.set_xlabel(xlabel)
    ax.set_xticks(np.arange(xticks[0], xticks[-1]+1))
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    if legend is True:
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.grid(True)
    return


def cst(rundir=None, damp=None, labels='auto', model=None, map=False,
        plot=True, it_no=10, verbose=False, plot_mzero=False, R=-0.2,
        long_title=True, SF=None, **plotargs):
    """
    rundir: path to rundir ('PATH/RUN'), motherdirectory of 'inversion_out'
            or
            list of paths
            or
            [PATH, (number of inversion_out folder)]
            or
            [
                [PATH1, (number of inversion_out folder)],
                [PATH2, (number of inversion_out folder)]
            ]

            or
                If SF is given, rundir=None

    damp: corresponding damping value
          or
          list of damping values
          or
          'all'

    it_no: number of iteration of the file
           or
           'all'
    """
    # If no splitting function is given, load data
    if SF is None:
        SF = read_cst_files(rundir, damp, labels, model,
                            it_no, verbose, plot_mzero, R,
                            long_title)
    if type(SF) is list or SplittingFunc:
        SF = Set(SF)

    if verbose is True:
        print(SF)

    if map is False and plot is True:
        SF.plot(R=R, **plotargs)
    elif plot is True:
        SF.plot_map(kind='cst', R=R, **plotargs)
        try:
            SF.plot_map(kind='dst', R=R, **plotargs)
        except Exception:
            pass

    return SF


def plot_coeffs_per_mode(ifiles=None, label1=None, SF_in=None, color1='r',
                         mode='sc', degree=2, order='all', l='all', n='all',
                         model=None, ifiles2=None, label2=None, color2='b',
                         plot_f_center=False, plot_Q=False, kind='cst',
                         fig_size=None, rot=90, ylim='auto', cmap='Dark2',
                         ordering=['n', 'l'], legend_show=True, fig_abc=None,
                         marker='o', markersize=8, spacing=None, R=-0.2,
                         yaxs_split=False, fs=10, savefig=False, verbose=False):
    """
    param ifiles: input files of first data set

    param label1: label of input files of first data set

    param color1: color of first data set

    param kind: cst or dst plot

    param mode: kind of modes to be plotted
    type  mode: string
    input mode: 'sc' - plots self-coupled coefficients
                'cc' - plots cross-coupled coefficients
                's'  - plots spheroidal mode coefficients
                't'  - plots toroidal mode coefficients

    param degree: s degree you want to plot. zero for fc and Qc plot

    param order: t degree you want to plot.

    param l: filter to plot modes with l=l

    param n: filter to plot modes with n=n

    param model: model predictions plotted with the data
    type  model: list of strings
    input model: available models:
                        'S20RTS', 'RR', 'AD', 'CB', 'DE', 'HT', 'QM1',
                        'REM', 'SAS', 'TCB', 'TZ'
                 e.g. ['S20RTS', 'RR']

    param ifiles2: input files of second data set to be pltted

    param label2: label of 2nd data set to be plotted

    param color2: color of second data set

    param fig_size: fig size

    param rot: rotate modes label

    param ylim: ylim in plot

    param cmap: cmap used to plot data

    param plot_f_center: if set 'True', degree 0 will be replaced with
                         the center frequency
    type  plot_f_center: bool
    input plot_f_center: True, False

    param plot_Q: plot Qc plot

    param ordering: specifies ordering of the modes in the branch, e.g.:
                    ['n', 'l'] orders first by n and then by s
    type  ordering: string, list
    input ordering: 'n', 'l', 'sens' or list of them e.g.
                    ['n', 'l'], sensitive to order

    param marker: set marker for data

    param markersize: size of marker to be plotted

    param spacing: spacing between markers plotted

    param ordering: order by type, n, l the branch plot

    param fig_abc: set the a), b) ... z) in the Fig title

    param fs: figure size

    param yaxs_split: split plot in yaxis to plot outlier separetely
                      if not False is should be a list:
                      [ymin_plot_bottom, ymax_plot_bottom,
                       ymin_plot_bottom, ymax_plot_top]

    param savefig: savefig

    param verbose: verbose

    param SF_in: labels should be in the same order as the SF in SF_in

    param R: To change r factor in dst prediction


    Example:
    Plots coefficients for a branch of modes, in this example
    all crosscoupling

    mdir = '/net/home/simons/mounts/eejit/splitting/modes/'
    mdir += '*-*/cross_coupling/'
    mdir += 'TZ-Comp/deg12/inversion_out/mnew-it9-d%4.2e.dat' % damp
    ifiles = glob.glob(mdir)
    ifiles.sort()

    for m in ifiles:
       if '00t19' in m:
           print(m)
           ifiles.remove(m)

    SF = plot_coeffs_per_mode(ifiles, model=['S20RTS', 'RR', 'CB'],
                              mode='sc', model=["HT", "QM1", "REM"],
                              plot_f_center=True,
                              degree = 0, ordering=["sens", "l"], l=2,
                              cmap='grey', rot=0,)

    """

    mpl.rcParams.update({'font.size': fs})
    if mode not in ['sc', 'cc', 's', 't']:
        if type(mode) != list:
            mode = [mode]
        mode = [format_name(x).upper() for x in mode]

    if plot_f_center is True:
        kind = 'cst'
    elif plot_Q is True:
        kind = 'dst'

    # Load Data here
    if SF_in is None:
        SF = Set()

        Nfiles = len(ifiles)
        if model is not None:
            Nmodels = Nfiles * len(model)
        else:
            model = [None]

        if ifiles2 is not None:
            if label2 is None:
                label2 = 'Data 2'
            model.append(label2)

        for i, f in enumerate(ifiles):
            if verbose is True:
                print('Loading Data: %s' % f)
            setup = join(f.split('inversion_out')[0], 'setup.pickle')
            setup = read_setup(setup)
            name = setup.pbsjob

            SF += loadmodel(setup=setup, ifile=f, format=None, damp=1,
                            name=name)
            if verbose is False:
                update_progress((i+1)/float(Nfiles), 'Loading   Data')

        if ifiles2 is not None:
            for f in ifiles2:
                if verbose:
                    print(f)
                setup = join(f.split('inversion_out')[0], 'setup.pickle')
                setup = read_setup(setup)
                name = setup.pbsjob
                SF += loadmodel(setup=setup, ifile=f, format=None, damp=1,
                                name="%s_%s" % (label2, name))

        for i, f in enumerate(ifiles):
            # counter += 1
            if verbose is True:
                print('Calculating Model for: %s' % f)
            setup = join(f.split('inversion_out')[0], 'setup.pickle')
            setup = read_setup(setup)
            name = setup.pbsjob

            for j, m in enumerate(model):
                if verbose is True:
                    print("Current model: %s" % m)
                if m is None or m == label2:
                    continue
                # counter += 1
                if verbose is False:
                    update_progress((i+1)/float(Nmodels), 'Loading Models')
                SF += loadmodel(setup=setup, ifile=None, format=m, damp=1,
                                name="%s_%s" % (m, name))

    if SF_in is not None:  # prevent input to be overwritten by output
        SF = SF_in.copy()

    if SF_in is not None and model is not None:
        # creating setup based on input SF
        for sf in SF_in:
            input = OrderedDict()

            # modes sc
            input['modes'] = OrderedDict()
            for m in sf.stats['modes_in']:
                smax = 0
                mname = m.name

                if mname in sf.cst:
                    for s, c in sf.cst[mname].iteritems():
                        if smax < int(s):
                            smax = int(s)
                    input['modes'].update([(mname, smax)])

            if kind == "dst" and not plot_Q:
                input['modes_sc_dst'] = OrderedDict()
                for m in sf.stats['modes_in']:
                    smax = 2
                    mname = m.name

                    if mname in sf.cst:
                        for s, c in sf.cst[mname].iteritems():
                            if smax < int(s):
                                smax = int(s)
                        input['modes_sc_dst'].update([(mname, smax)])

            if len(sf.stats['modes_cc_in']) > 0:
                input['modes_cc'] = OrderedDict()
                for m in sf.stats['modes_cc_in']:
                    mname = m.name
                    cc = []
                    for m1 in mname.split('-'):
                        cc.extend(re.split(r'(S+|s+|T+|t+)', m1))
                    smin = min(max_cc_degrees(cc))
                    smax = 0
                    if mname in sf.cst:
                        for s, c in sf.cst[mname].iteritems():
                            if smax < int(s):
                                smax = int(s)
                        input['modes_cc'].update([(mname, [smin, smax])])
            setup = Setup('CST', input)
            counter = 0
            Nfiles = len(model)
            name = sf.stats.name
            for j, m in enumerate(model):
                counter += 1
                update_progress(counter/float(Nfiles), 'Loading Models')
                SF += loadmodel(setup=setup, ifile=None, format=m, damp=1,
                                name="%s_%s" % (m, name), R=R)

    if model is None:
        model = [None]
    # reading second set of SF types. Only two possible tight now
    m = []
    if verbose is True:
        print('\nLoaded Data/models:')
    for S in SF:
        if verbose is True:
            print(S)

        if S.stats.name.split('_')[0] not in m:
            m.append(S.stats.name.split('_')[0])

    if SF_in is not None:
        if model[0] is not None: # I dont know why this is here
            model.append(m[1])
        if len(m) == 2 and model[0] is None:
            model = [m[-1]] # apppend 2nd data set as 1st model

    # spacing between coeffs for the same modes,
    # only if one than one data set is plotted
    if spacing and model[0] is not None:
        input = list(model)
        input.insert(len(input), 'data')
        n_input = len(input)
        ww = np.linspace(-0.25, 0.25, n_input)
        width = {}
        for m, w in zip(input, ww):
            width[m] = w
    else:
        input = list(model)
        input.insert(len(input), 'data')
        width = {}
        for m in input:
            width[m] = 0

    # Prepare Figures and Plot
    # use first splitting function to check the amount of plots
    # we need 2 * degree + 1 figures
    N = 2 * degree + 1
    rows = int(np.ceil(N/2.0))
    if order is 'all' and plot_Q is False and plot_f_center is False:
        # cst: plot all t orders of a certain degree
        cols = 2
        if fig_size is None:
            fig_size = (7.5, round(0.75*N))
        gs = gridspec.GridSpec(rows, cols,  wspace=0.1)
        fig = plt.figure()
        fig.set_size_inches(fig_size)

    elif order is 0 or plot_Q is True or plot_f_center is True:
        cols = 1
        if mode is 'sc' and fig_size is None:
            fig_size = (7, 2)
        elif mode is 'cc' and fig_size is None:
            fig_size = (7, 2.5)
        gs = gridspec.GridSpec(1, 1)
        fig = plt.figure()
        fig.set_size_inches(fig_size)

    elif plot_Q is True or plot_f_center is True:
        cols = 1
        gs = gridspec.GridSpec(1, 2)
        fig = plt.figure()
        fig.set_size_inches(fig_size)

    # Prepare markers and ticks
    ax_dict = {}
    s_nums = {}
    label_model_set = {}
    mark = {}
    x_ticklabel = []
    i = 0

    nosym = [",", ".", "<", "X", "x", "P", "D", "+", "*",
             "tri_down", "tri_up", "tri_left", "tri_right", "plus",
             "vline", "hline", "tickdown", "tickup", "tickleft", "tickright",
             "1","2","3","4","p","v"]
    markers = mmarkers.MarkerStyle().markers
    markers = {k: v for k, v in markers.items() if v != 'nothing'}
    markers = {k: v for k, v in markers.items() if type(k) is not int}
    # markers = {k: v for k, v in markers.items() if k.isdigit() is not True}
    markers = {k: v for k, v in markers.items() if k not in nosym}
    markers = iter(markers)
    colors = {}
    if cmap is not 'grey':
        if label2 is not None:
            colormap = get_iter_colormap(model[0:-1], cmap,
                                         random_values='center')
        else:
            colormap = get_iter_colormap(model, cmap, random_values='center')
    # Prepare xtick-labels and colors

    if verbose is True:
        print('\nPreparing labels for: ')
    for s in SF:
        if verbose is True:
            print(s)
        name = s.stats.name.split('_')[0]
        mname = None
        if name in model:
            mname = name
            try:
                name = s.stats.name.split('_')[1]
            except Exception:
                name = s.stats.name
        if mode == 'cc':
            labels = s.stats.modes_cc_in.names
            if verbose is True:
                msg = 'Reading labels: %s' % labels
                print(msg)
            m_list = labels[:]  # removing modes from branch
            for m in m_list:
                mm = split_digit_nondigit(m)
                if l is not 'all' and isinstance(l, (list,)):
                    if int(mm[2]) is l[0] and int(mm[6]) is l[1]:
                        pass
                    elif int(mm[2]) is l[1] and int(mm[6]) is l[0]:
                        pass
                    else:
                        labels.remove(m)
            labels = [format_name(x, 2) for x in labels]
        else:
            labels = s.stats.modes_in.names
            if verbose is True:
                msg = 'Reading labels: %s' % labels
                print(msg)
            m_list = labels[:]  # removing modes from branch
            for m in m_list:
                mm = split_digit_nondigit(m)
                if l is not 'all' and isinstance(l, (int,)):
                    if int(mm[2]) is not l:
                        labels.remove(m)
                if n is not 'all':
                    if int(mm[0]) is not n:
                        labels.remove(m)
            labels = [format_name(x, 2) for x in labels]

        for label in labels:
            if label not in x_ticklabel and (mname is None or mname==label2):
                if mode in ['s', 't', 'S', 'T']:
                    if mode.upper() not in label:
                        continue
                x_ticklabel.append(label)
                s_nums[label] = i
                i += 1
            else:
                mlabel = mname
                label_model_set[mlabel] = False
                if mlabel == label2:
                    mark[mlabel] = marker
                    colors[mlabel] = color2
                else:
                    if mlabel not in mark:
                        #if mname=="S20RTS":
                        #    mark[mlabel] = 'D'
                        #else:
                        #    mark[mlabel] = 's'#markers.next()  #
                        mark[mlabel] = markers.next()  #
                    if mlabel not in colors:
                        if cmap == 'grey':
                            colors[mlabel] = 'grey'
                        else:
                            try:
                                colors[mlabel] = colormap.next()
                            except StopIteration:
                                if label2 is not None:
                                    colormap = get_iter_colormap(
                                                        model[0:-1], cmap,
                                                        random_values='center')
                                else:
                                    colormap = get_iter_colormap(
                                                        model, cmap,
                                                        random_values='center')
                                colors[mlabel] = colormap.next()

    # Ordering modes
    mode_list = Modes()
    labels = []
    if mode != 'cc':  # SC splittingfunc
        for x in x_ticklabel:
            if verbose:
                msg = 'Reading mode for tick: %s' % x
                print(msg)
            m = read_modes(modenames=str(x))[0]
            mode_list += m
        if isinstance(ordering, (list,)):
            mode_list.sort(keys=ordering)
        else:
            mode_list.sort(keys=[ordering])
        for i, tick in enumerate(mode_list):
            if verbose:
                msg = 'Setting tick: %s' % tick
                print(msg)
            x_ticklabel[i] = "$_{%s}%s_{%s}$" % (tick.n, tick.type, tick.l)
            labels += [format_name(tick.name, 2)]
            s_nums[format_name(tick.name, 2)] = i

    else:  # CC splittingfunc
        for x in x_ticklabel:
            m = Mode()
            m.name = x

            # CC sens the same as 2nd mode in CC
            x = x.split('-')[1]
            m_sc = read_modes(modenames=str(x))[0]
            m.sens = m_sc.sens
            mode_list += m

        if isinstance(ordering, (list,)):  # only works for sens
            mode_list.sort(keys=ordering)
        else:
            mode_list.sort(keys=["cc"])

        for i, tick in enumerate(mode_list):
            m = split_digit_nondigit(tick.name)
            x_ticklabel[i] = ("$_{%s}%s_{%s}$-$_{%s}%s_{%s}$" %
                              (int(m[0]), m[1], int(m[2]),
                               int(m[4]), m[5], int(m[6])))
            labels += [format_name(tick.name, 2)]
            s_nums[format_name(tick.name, 2)] = i

    # Prepare Vlines that separates overtones here
    x_init = range(len(s_nums))
    label_data_set = False
    # Find ticks where put a vline
    # e.g. 00T04 and 01t04, should have a vline, to seperate the branches
    ii = 1
    set_line = []

    if verbose is True:
        print('\nPreparing vlines for:')
    for x, y, m1, m2 in zip(labels[:], labels[1:],
                            mode_list[:], mode_list[1:]):

        if mode != 'cc':  # SC splittingfunc
            if "n" in ordering:
                n1 = int(split_digit_nondigit(x)[0])
                n2 = int(split_digit_nondigit(y)[0])
            elif "l" in ordering:
                n1 = int(split_digit_nondigit(x)[2])
                n2 = int(split_digit_nondigit(y)[2])
            if "sens" in ordering:
                n1 = m1.sens
                n2 = m2.sens
            if "freq" in ordering:
                n1 = int(m1.freq)
                n2 = int(m2.freq)
            if n1 != n2:
                set_line.append(ii)
                if verbose is True:
                    print("Line between %s %s" % (x, y))
            ii += 1
        else:  # CC splittingfunc
            if "n" in ordering:
                n1 = int(split_digit_nondigit(x)[0])
                n2 = int(split_digit_nondigit(y)[0])
                n3 = int(split_digit_nondigit(x)[4])
                n4 = int(split_digit_nondigit(y)[4])
            elif "l" in ordering:
                n1 = int(split_digit_nondigit(x)[2])
                n2 = int(split_digit_nondigit(y)[2])
                n3 = int(split_digit_nondigit(x)[6])
                n4 = int(split_digit_nondigit(y)[6])
            if "sens" in ordering:
                n1 = m1.sens
                n2 = m2.sens
                n3 = m1.sens
                n4 = m2.sens
            if n1 != n2 or n3 != n4:
                set_line.append(ii)
                if verbose is True:
                    print("Line between %s %s" % (x, y))
            ii += 1

    # Plotting loop
    if verbose is True:
        print('\nPlotting:')
    for s in SF:
        if kind is "cst":
            # cst
            sf = s.cst
            sf_err = s.cst_errors
        elif kind is 'dst':
            # dst
            sf = s.dst
            sf_err = s.dst_errors

        for key in sf.iterkeys():
            skip = False
            if mode not in ['cc', 'sc']:
                if '-' in key:
                    continue
                for m in mode:
                    if m.upper() not in format_name(key):
                        skip = True
            elif mode == 'cc' and '-' not in key:
                continue
            elif mode == 'sc' and '-' in key:
                continue
            # Only plot selected modes
            if mode == 'sc':
                m = split_digit_nondigit(key)
                if l is not 'all' and isinstance(l, (int,)):
                    if int(m[2]) is not l:
                        skip = True
                if n is not 'all':
                    if int(m[0]) is not n:
                        skip = True

            if mode == 'cc':
                m = split_digit_nondigit(key)
                if l is not 'all' and isinstance(l, (list,)):
                    if int(m[2]) is l[0] and int(m[6]) is l[1]:
                        pass
                    elif int(m[2]) is l[1] and int(m[6]) is l[0]:
                        pass
                    else:
                        skip = True
            if skip is True:
                continue
            if str(degree) not in sf[key] and degree not in sf[key]:
                continue

            if order is 'all':
                coeffs = sf[key][str(degree)]
                try:
                    errors_temp = sf_err[key][str(degree)]
                    errors = errors_temp['uncertainty']
                    errors_up = errors_temp['upper_uncertainty']
                    errors_lw = errors_temp['lower_uncertainty']
                    if errors_up is None:
                        errors_up = [None]*len(coeffs)
                        errors_lw = [None]*len(coeffs)
                except Exception:
                    errors = [None]*len(coeffs)
                    errors_up = [None]*len(coeffs)
                    errors_lw = [None]*len(coeffs)

            elif order is 0:
                coeffs = [sf[key][str(degree)][0]]
                try:
                    errors_temp = sf_err[key][str(degree)]
                    errors = [errors_temp['uncertainty'][0]]
                    errors_up = [errors_temp['upper_uncertainty'][0]]
                    errors_lw = [errors_temp['lower_uncertainty'][0]]
                except Exception:
                    errors = [None]
                    errors_up = [None]
                    errors_lw = [None]

            for i, (cst, err, erru, errl) in enumerate(zip(coeffs, errors,
                                                           errors_up,
                                                           errors_lw)):
                eval = err
                # Check if degree is 0, then plot fc/Q
                if (degree == 0 and (plot_f_center is True or plot_Q is True)):
                    cst_fQ = None
                    if verbose:
                        print(key, type(key))
                    _mode = read_modes(modenames=str(key))[0]
                    if plot_f_center is True:
                        # fc = f0 + (4pi)**-1/2 * Re(c00) in microHerz
                        c00 = s.cst[key]['0']
                        fc  = _mode.freq * 1e3 + 1./np.sqrt(4. * np.pi) * c00
                        cst = (fc - _mode.freq * 1e3)
                        if verbose is True and err is None:
                            print(r'c00=%s, f=%s' % (c00, cst))
                        if err is not None:
                            if verbose is True and erru is None:
                                print('c00=%.2f +/- %.2f' % (c00, err))
                            err = 1./np.sqrt(4. * np.pi) * err
                            if verbose is True and erru is None:
                                print('f=%.2f +/- %.2f' % (fc, err))
                        if erru is not None:
                            if verbose is True:
                                print('c00=%.2f %.2f/ %.2f' % (c00, errl, erru))
                            erru = 1./np.sqrt(4. * np.pi) * abs(erru)
                            errl = 1./np.sqrt(4. * np.pi) * abs(errl)
                            if verbose is True:
                                print('f=%.2f -%.2f/ +%.2f\n' % (fc,
                                                                 errl, erru))
                            eval = np.array([errl, erru]).reshape((2, 1))
                        if plot_Q is True:
                            cst_fQ = {'f': cst}

                    if plot_Q is True:
                        # Q = fc / (2 * (f0/(2*Q0)+(4pi)**-1/2 * Im(c00)))
                        c00 = s.cst[key]['0']
                        fc = _mode.freq * 1e3 + 1. / np.sqrt(4. * np.pi) * c00
                        # fcN, QcN = cst2fQ(s.cst[key]['0'], s.dst[key]['0'],
                        #                   _mode)
                        Qc = calc_Q(_mode, fc, cst)

                        if verbose is True and err is None:
                            print('c00=%s, f=%s' % (c00[0], fc[0]))
                            print('d00=%s, Q=%s' % (cst, Qc[0]))

                        if err is not None:
                            if verbose is True and erru is None:
                                print('c00=%.2f, df=%.2f' % (c00[0], fc[0]))
                                print('d00=%.2f +/- %.2f' % (cst, err))
                            Qe = calc_Q(_mode, fc, cst, err)
                            err = Qc - Qe
                            if verbose is True and erru is None:
                                print('Q=%.2f +/- %.2f' % (Qc[0], err[0]))

                        if erru is not None:
                            # err in d00 translates to opposite sign error in Q
                            if verbose is True:
                                print('c00=%.2f, df=%.2f' % (c00[0], fc[0]))
                                print('d00=%.2f %.2f/ %.2f' % (cst, errl, erru))

                            if errl > 0:
                                # if its positive only upper uncsrt exist
                                errl = 0

                            if erru < 0:
                                # if its negative only lower uncsrt exist
                                erru = 0

                            Qe = calc_Q(_mode, fc, cst, erru)
                            erru = abs(Qc - Qe)
                            Qe = calc_Q(_mode, fc, cst, errl)
                            errl = abs(Qc - Qe)

                            # they are flipped on purpose
                            if verbose is True:
                                print('Q=%.2f -%.2f/ +%.2f\n' % (Qc[0],
                                                                 erru, errl))
                            eval = np.array([erru, errl]).reshape((2, 1))
                        cst = Qc - _mode.Q
                        if plot_f_center is True:
                            cst_fQ = {'Q': cst}

                # Find correct plot row and column here
                if i > 0:
                    j = i + 1  # skip subplot (0,1)
                else:
                    j = i
                i_row = int(np.ceil((j + 1)/2. - 1))
                if i_row < 0:
                    i_row = 0
                i_col = int(abs(np.sin(j*np.pi/2.)))

                if i not in ax_dict:
                    if i == 0:
                        if plot_f_center is True:
                            title = 'Center frequency'
                            weight="bold"
                        elif plot_Q is True:
                            title = 'Quality factor'
                            weight="bold"
                        elif kind is 'cst':
                            title = '$c_{%s0}$ coefficients of IC sensitive cross-coupling modes' % (degree)
                            weight=None
                        elif kind is 'dst':
                            title = '$d_{%s0}$' % (degree)
                            weight=None
                    else:
                        t = int(np.ceil(i/2.))
                        if i_col == 0:
                            if kind is 'cst':
                                title = '$\mathrm{Re}(c_{%s%s})$' % (degree, t)
                                weight=None
                            elif kind is 'dst':
                                title = '$\mathrm{Re}(d_{%s%s})$' % (degree, t)
                                weight=None
                        else:
                            if kind is 'cst':
                                title = '$\mathrm{Im}(c_{%s%s})$' % (degree, t)
                                weight=None
                            elif kind is 'dst':
                                title = '$\mathrm{Im}(d_{%s%s})$' % (degree, t)
                                weight=None

                    if isinstance(yaxs_split, list):
                        gridspec_kw = {'height_ratios': [0.3,1]}
                        fig, (axt2, axt) = plt.subplots(2,1,sharex=True,
                                                      gridspec_kw=gridspec_kw,
                                                      facecolor='w')
                        fig.set_size_inches(fig_size)
                    else:
                        axt = plt.subplot(gs[i_row, i_col])

                    if fig_abc:
                        if isinstance(yaxs_split, list):
                            axt2.set_title('%s) %s' % (fig_abc,title),
                                          weight=weight)
                        else:
                            axt.set_title('%s) %s' % (fig_abc,title),
                                          weight=weight)
                    else:
                        if isinstance(yaxs_split, list):
                            axt2.set_title(title, weight="bold")
                        else:
                            axt.set_title(title, weight="bold")

                    if plot_f_center is True or plot_Q is False:
                        if isinstance(yaxs_split, list):
                            axt2.set_ylabel('$f$ ($\mu$Hz)')
                        else:
                            axt.set_ylabel('$f$ ($\mu$Hz)')
                    if plot_Q is True:
                        if isinstance(yaxs_split, list):
                            axt2.set_ylabel('$\delta Q$')
                        else:
                            axt.set_ylabel('$\delta Q$')

                    # axt.grid()
                    axt.hlines(0, -1, len(s_nums), color='lightgray')
                    axt.set_xlim(-0.5, len(s_nums)-0.5)
                    axt.set_xticks(x_init)
                    if i_row == rows-1 or order is 0:
                        axt.set_xticklabels(x_ticklabel, fontsize=12)
                        for tick in axt.get_xticklabels():
                            tick.set_rotation(rot)
                    else:
                        axt.tick_params(labelbottom=False)

                    if ylim != 'auto':
                        axt.set_ylim(ylim)

                    # split y axis
                    if isinstance(yaxs_split, list):
                        axt.set_ylim(yaxs_split[0], yaxs_split[1])
                        axt2.set_ylim(yaxs_split[2], yaxs_split[3])

                        axt.spines['top'].set_visible(False)
                        axt2.spines['bottom'].set_visible(False)
                        axt.tick_params(top=False)
                        axt2.tick_params(bottom=False)

                        # directly from:
                        # https://matplotlib.org/examples/pylab_examples/broken_axis.html
                        d = .01  # how big to make the diagonal lines in ax
                        kwargs = dict(transform=axt.transAxes, color='k',
                                      clip_on=False)
                        axt.plot((-d,+d),(1-d,1+d),**kwargs) # bottom/left
                        axt.plot((1-d,1+d),(1-d,1+d),**kwargs) # bottom/right

                        kwargs.update(transform=axt2.transAxes, color='k',
                                      clip_on=False)
                        axt2.plot((1-d,1+d),(-d,+d),**kwargs) # top/right
                        axt2.plot((-d,+d),(-d,+d),**kwargs) # top/left
                    ax_dict[i] = axt

                ax = ax_dict[i]
                name = format_name(key, 2)
                mname = s.stats.name.split('_')[0]

                if err is None:
                    eval = None
                if mname in model:
                    if verbose is True:
                        print('current model: %s | %s' % (mname, name))
                    try:
                        s_num = s_nums[name]
                    except KeyError:
                        continue

                    zorder = 10
                    if mname == label2:
                        zorder = 150
                    else:
                        eval = err

                    if label_model_set[mname] is False:
                        ax.errorbar(s_num+width[mname], cst, yerr=eval,
                                    marker=mark[mname],
                                    markersize=markersize, color=colors[mname],
                                    linestyle='None', elinewidth=2,
                                    capsize=3, zorder=zorder, label=mname)
                        label_model_set[mname] = True
                        # removing whiskers from legend
                        ehandles, elabels = ax.get_legend_handles_labels()
                        ehandles = [h[0] for h in ehandles]
                        if verbose is True:
                            msg = 'current model: %s | %s - plotted, label'
                            msg = msg % (mname, name)
                            print(msg)
                        if N != 1 and order is 'all':
                            x_pos = cols - 0.9
                            y_pos = 1
                            if legend_show:
                                ax.legend(ehandles, elabels, frameon=False,
                                        handletextpad=0.25,
                                        bbox_to_anchor=(x_pos, y_pos),
                                        loc=2, borderaxespad=0.)
                        else:
                            if legend_show:
                                ax.legend(ehandles, elabels, frameon=False,
                                        handletextpad=0.25,
                                        bbox_to_anchor=(1, 0.5),
                                        loc='center left', borderaxespad=0.)
                    else:
                        ax.errorbar(s_num+width[mname], cst, yerr=eval,
                                    marker=mark[mname],
                                    markersize=markersize, color=colors[mname],
                                    linestyle='None', elinewidth=2,
                                    capsize=3, zorder=zorder)
                        if isinstance(yaxs_split, list):
                            axt2.errorbar(s_num+width[mname], cst, yerr=eval,
                                        marker=mark[mname],
                                        markersize=markersize, color=colors[mname],
                                        linestyle='None', elinewidth=2,
                                        capsize=3, zorder=zorder)
                        if verbose is True:
                            msg = 'current model: %s | %s - plotted'
                            msg = msg % (mname, name)
                            print(msg)
                else:
                    if verbose is True:
                        print('current data: %s' % name)
                    s_num = s_nums[name]
                    # Check for current n of mode, if new branch starts
                    # draw a vertical line
                    if s_num in set_line:
                        ax.axvline(x=s_num-0.5, color='lightgray',
                                   linestyle='--')
                        if isinstance(yaxs_split, list):
                            axt2.axvline(x=s_num-0.5, color='lightgray',
                                       linestyle='--')
                    if label_data_set is False:
                        if label1 is None:
                            label1 = 'Data'
                        ax.errorbar(s_num+width['data'], cst, yerr=eval,
                                    marker=marker,
                                    markersize=markersize, color=color1,
                                    linestyle='None', elinewidth=2,
                                    capsize=3, zorder=100, label=label1,)
                        label_data_set = True
                        # removing whiskers from legend
                        ehandles, elabels = ax.get_legend_handles_labels()
                        ehandles = [h[0] for h in ehandles]

                        if verbose is True:
                            print('current data: %s - plotted, label' % name)

                        if N != 1 and order is 'all':
                            x_pos = cols - 0.9
                            y_pos = 1
                            if legend_show:
                                ax.legend(ehandles, elabels, frameon=False,
                                        handletextpad=0.25,
                                        bbox_to_anchor=(x_pos, y_pos),
                                        loc=2, borderaxespad=0.)
                            else:
                                if legend_show:
                                    ax.legend(ehandles, elabels, frameon=False,
                                        handletextpad=0.25,
                                        bbox_to_anchor=(1, 0.5),
                                        loc='center left', borderaxespad=0.)
                    else:
                        ax.errorbar(s_num+width['data'], cst, yerr=eval,
                                    marker=marker,
                                    markersize=markersize, color=color1,
                                    linestyle='None', elinewidth=2,
                                    capsize=3, zorder=100)
                        if verbose is True:
                            print('current data: %s - plotted' % name)
    if savefig:
        if degree is 0 or order is 0:
            t = 0
        elif order is 'all':
            t = 't'
        fname = '%s_%s%s%s' % (mode, kind[0], degree, t)
        plt.tight_layout()
        fig.savefig('cbranch_%s.png' % fname, orientation='landscape', dpi=400,
                    bbox_inches='tight', pad_inches=0.025, )#transparent=True)

    return SF


def hsnr(schema, damping, setup=None, run_dir=None, iteration=10,
         verbose=False, fig_size='auto'):
    """
    Added hsnr according to Resovsky and Ritzwoller 1998
    """
    if setup is None and run_dir is None:
        raise IOError('No setup or run_dir provided')
    if setup is not None:
        run_dir = join(setup.rundir, setup.inversion_outdir)
    if setup is None and run_dir is not None:
        setup = read_setup(join(run_dir, 'setup.pickle'))
    db_path = get_db_path(run_dir)
    if db_path is None:
        return

    fig, ax = plt.subplots()
    for i, s in enumerate(schema):
        # Query Inv(ersion) object
        hsnr = []
        Inv = query_ev_seg_misfits(db_path, s, [iteration], [damping],
                                   per_station=True)
        events = Inv.allmisfits[damping].keys()
        seg = read_seg('db', channel=s, modes=setup.modes.names, events=events)

        for cmt, values in Inv.allmisfits[damping].iteritems():
            for station, mf in values:
                try:
                    snr = seg.select(station=station, event=cmt)[0].stats.snr
                    hsnr.append([mf, snr])
                except Exception:
                    if verbose:
                        print(s, cmt, station)
                    continue
        hsnr = np.array(hsnr)
        ax.scatter(hsnr.transpose()[1], hsnr.transpose()[0], 5, label=s)
    ax.set_ylabel('Misfit')
    ax.set_xlabel('SNR')
    if fig_size != 'auto':
        fig.set_size_inches(fig_size)
    plt.legend()

    return


def calc_Q(mode, fc, cst, err=None):
    if err is None:
        _Q = ((mode.freq * 1e3) / (2*mode.Q)) + \
              1. / np.sqrt(4. * np.pi) * cst
    else:
        _Q = ((mode.freq * 1e3) / (2*mode.Q)) + \
              1. / np.sqrt(4. * np.pi) * (cst + err)
    return 0.5 * fc / _Q
