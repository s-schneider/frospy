import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.markers as mmarkers
from matplotlib.pyplot import cm

import re
from collections import OrderedDict

from frospy.splitting.load import load
from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes, Mode, format_name, calc_Q
from frospy.core.setup.settings import Setup
from frospy.util.base import (update_progress, split_digit_nondigit,
                              max_cc_degrees)
from frospy.plot.nmplt import get_iter_colormap, multicolor_ylabel


def branch(ifiles=None, data_label=None, label1=None, SF_in=None,
           color1='red',
           mode='sc', degree=2, order='all', l='all', n='all',
           model=None, ifiles2=None, label2=None, color2='b',
           plot_f_center=False, plot_Q=False, kind='cst',
           fig_size=(12, 5), rot=90, ylim='auto', cmap='Dark2',
           cmap_random_values='center',
           model_colors=None,
           ordering=['n', 'l'], legend_show=True, fig_abc=None,
           marker='o', markersize=8, spacing=None,
           labeltop=False, R=-0.2,
           yaxs_split=False, fs=22, savefig=False,
           border_width=3, tick_width=1, format_yticks=None,
           errorbars='asymmetrical', linestyle=None,
           elinewidth=2, wspace=0.1,
           verbose=False, loadingbar=False, damping_label=None,
           bbox_inches='tight', pad_inches=0.025,
           cmap_all_damping='Reds_r', suptitle=None,
           linewidth_S20=2, add_colorbar=False,
           cbar_aspect=1.2,
           vmin=-4, vmax=4,
           colorbarlabel=None,
           include_CRUST=True,
           colorbar_multicolor=['blue', 'k', 'red'],
           colorbarlabel_anchor=(0.5, -.8),
           subtitle=None,
           dpi=400,
           exclude=None,
           parts='all',
           skip2=-True,
           **savefigargs):
    """
    ATTENTION: This is the WIP version of plot_coeffs_per_mode,
    If you are unsure, please use plot_coeffs_per_mode

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

    param include_CRUST: add CRUST5.1 to model predictions
    type  include_CRUST: bool (add to all models) or
                         list of models where it will be added
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

    param border_width: sets the line width of the frame of the plot

    param labeltop: bool, if True modelabels are printed on the top side too.

    param fs: figure size

    param yaxs_split: split plot in yaxis to plot outlier separetely
                      if not False is should be a list:
                      [ymin_plot_bottom, ymax_plot_bottom,
                       ymin_plot_bottom, ymax_plot_top]

    param savefig: savefig

    param verbose: verbose

    param SF_in: labels should be in the same order as the SF in SF_in

    param R: To change r factor in dst prediction

    param errorbars: 'symmetrical' or 'asymmetrical'

    Example:
    Plots coefficients for a branch of modes, in this example
    all crosscoupling

    mdir = '/tmp/eejit_simons/splitting/modes/'
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

    if mode not in ['sc', 'cc', 's', 't', 'stoneley']:
        if type(mode) != list:
            mode = [mode]
        mode = [format_name(x).upper() for x in mode]

    if plot_f_center is True:
        kind = 'cst'
        degree = 0
    elif plot_Q is True:
        kind = 'dst'
        degree = 0

    if data_label is None:
        data_label = ['Data']

    # check if data is in SF
    if SF_in is not None:  # prevent input to be overwritten by output
        SF = SF_in.copy()

    # If no data is found, the first model is set to data
    data_found = False
    for S in SF:
        if S.stats.name == 'data':
            data_found = True
    if data_found is not True:
        ref_name = SF[0].stats.name
        data_label = [ref_name]
        for _i, S in enumerate(SF):
            if S.stats.name == ref_name:
                SF[_i].stats.name = 'data'

    if SF_in is not None and model is not None:
        # creating setup based on input SF
        counter = 0
        Nfiles = len(SF_in) * len(model)
        for sf in SF_in:
            input = OrderedDict()

            # modes sc
            input['modes'] = OrderedDict()
            for m in sf.stats['modes_in']:
                smax = 0
                mname = m.name

                if mname in sf.cst:
                    for s, c in sf.cst[mname].items():
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

            if len(sf.stats['modes_cc_in']) > 0 and sf.stats['modes_cc_in'] != 'None':
                input['modes_cc'] = OrderedDict()
                for m in sf.stats['modes_cc_in']:
                    mname = m.name
                    cc = []
                    for m1 in mname.split('-'):
                        cc.extend(re.split(r'(S+|s+|T+|t+)', m1))
                    smin = min(max_cc_degrees(cc))
                    smax = 0
                    if mname in sf.cst:
                        for s, c in sf.cst[mname].items():
                            if smax < int(s):
                                smax = int(s)
                        input['modes_cc'].update([(mname, [smin, smax])])
            setup = Setup('CST', input)
            if SF_in is None:
                counter = 0
                Nfiles = len(model)

            name = sf.stats.name

            if name.lower() != 'data' and name.lower().startswith('data'):
                continue

            for j, m in enumerate(model):
                if include_CRUST is True:
                    CRUST = True
                else:
                    if type(include_CRUST) == list and model in include_CRUST:
                        CRUST = True
                    else:
                        CRUST = False
                counter += 1
                if loadingbar is True:
                    update_progress(counter/float(Nfiles), 'Loading Models')
                else:
                    if verbose is True:
                        print('Calculating model / modes: %s / %s' % (m, sf))
                SF += load(setup=setup, ifile=None, format=m, damp=1,
                           name="%s_%s" % (m, name), R=R,
                           include_CRUST=CRUST)

    if model is None:
        model = [None]
    # reading second set of SF types. Only two possible tight now
    m = []
    datain = {}
    if verbose is True:
        print('\nLoaded Data/models:')
    for S in SF:
        if verbose is True:
            print(S)
        if S.stats.name != 'data' and S.stats.name.startswith('data'):
            name = S.stats.name.split()[0]
            if name not in datain:
                datain[name] = S.stats.model

        if S.stats.model not in m:
            if S.stats.name == 'data':
                # print(S.stats.name, S.stats.model)
                m.append(S.stats.model)
            elif not S.stats.name.startswith('data'):
                # print(S.stats.name, S.stats.model)
                m.append(S.stats.model)

    if verbose is True:
        print('models', m)
    if SF_in is not None:
        # if model[0] is not None: # I dont know why this is here
        #     model.append(m[1])
        if len(m) == 2 and model[0] is None:
            model = [m[-1]]  # apppend 2nd data set as 1st model
        else:
            model = m
    if verbose is True:
        print('models w/o data', model)
    # spacing between coeffs for the same modes,
    # only if one than one data set is plotted
    from IPython import embed; embed()
    dkey = 'data'
    if spacing: # and
        if model[0] is None:
            input = []
            for d in datain.keys():
                input.append(d)
            if 'data' in input:
                dkey = 'measurement'
            input.append(dkey)
            for _m in model[1:]:
                input.append(_m)

        else:
            input = list(model)

        n_input = len(input)
        # right now S20, S40 and SP12 are plotted as a line by Default
        # if we make that optional we have to change this if condition too
        # If S20/S40 are not lines, but dots, they need to be taken
        # into account, as well as 'data'
        for xin in ('S20RTS', 'S40RTS', 'SP12RTS'):
            if xin in input:
                n_input += -1

        if len(datain) == 0:
            n_input += -1

        if type(spacing) in (int, float):
            ww = np.linspace(-spacing, spacing, n_input)
        else:
            ww = np.linspace(-0.25, 0.25, n_input)

        width = {}
        i = 0
        for m in input:
            # right now S20 and S40 are plotted as a line by Default
            # if we make that optional we have to change this if condition too
            if m in ('S20RTS', 'S40RTS', 'SP12RTS'):
                width[m] = 0
            if m == dkey and len(datain) == 0:
                width[m] = 0
            else:
                width[m] = ww[i]
                i += 1
    else:
        input = list(model)
        input.insert(-1, dkey)
        width = {}
        for m in input:
            width[m] = 0

    if verbose is True:
        print('Final models ', model)
        print('Spacing ', width)
    # Prepare Figures and Plot
    # use first splitting function to check the amount of plots
    # we need 2 * degree + 1 figures
    if parts == 'all':
        N = 2 * degree + 1
    else:
        N = 1
    rows = int(np.ceil(N/2.0))
    if add_colorbar == 'bottom':
        rows += 1

    if order == 'all' and plot_Q is False and plot_f_center is False and parts == 'all':
        # cst: plot all t orders of a certain degree
        cols = 2
        if add_colorbar == 'side':
            cols += 1
        if fig_size is None:
            fig_size = (7.5, round(0.75*N))
        if add_colorbar is False:
            gs = gridspec.GridSpec(rows, cols,  wspace=wspace)
        elif add_colorbar == 'bottom':
            wr = (rows-1) * [1]
            wr += [0.08]
            gs = gridspec.GridSpec(rows, cols,  wspace=wspace,
                                   height_ratios=wr)
        elif add_colorbar == 'side':
            wr = (cols-1) * [1]
            wr += [0.05]
            gs = gridspec.GridSpec(rows, cols,  wspace=wspace,
                                   width_ratios=wr)
        fig = plt.figure()
        fig.set_size_inches(fig_size)

    elif order == 0 or plot_Q is True or plot_f_center is True or parts != 'all':
        cols = 1
        if mode == 'sc' and fig_size is None:
            fig_size = (7, 2)
        elif mode == 'cc' and fig_size is None:
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

    nosym = [",", ".", "<", "X", "x", "P", "D", "+", "o",
             "tri_down", "tri_up", "tri_left", "tri_right", "plus",
             "vline", "hline", "tickdown", "tickup", "tickleft", "tickright",
             '^', '>', "1", "2", "3", 'H']
    markers = mmarkers.MarkerStyle().markers
    markers = {k: v for k, v in markers.items() if v != 'nothing'}
    markers = {k: v for k, v in markers.items() if type(k) is not int}
    # markers = {k: v for k, v in markers.items() if k.isdigit() is not True}
    markers = {k: v for k, v in markers.items() if k not in nosym}

    marker_order = ['d', 's', 'h', 'v', '*', 'p', '8', '4']
    _markers = OrderedDict()
    for m in marker_order:
        _markers.update([(m, markers[m])])

    if verbose:
        print('models :', len(model), 'markers :', len(markers))
    markersit = iter(_markers)
    markers = {}
    for _m in model:
        if _m is None:
            continue
        try:
            markers[_m] = next(markersit)
        except StopIteration:
            markersit = iter(_markers)
            markers[_m] = next(markersit)
        if verbose is True:
            print('markers: ', _m, markers[_m])
    for _d in datain.values():
        try:
            markers[_d] = next(markersit)
        except StopIteration:
            markersit = iter(_markers)
            markers[_d] = next(markersit)
        if verbose is True:
            print('markers: ', _d, markers[_d])

    colors = {}
    colors_data = {}

    # _label_damping is set later, before the plotting loop. Maybe it can be
    # moved here too
    _label_damping = []
    if cmap != 'grey':
        if label2 is not None:
            colormap = get_iter_colormap(model[0:-1], cmap,
                                         random_values=cmap_random_values)
        else:
            colormap = get_iter_colormap(model, cmap,
                                         random_values=cmap_random_values)

    _data = []
    for s in SF:
        if s.stats.model is None:
            _data.append(s.stats.name)

    _data = list(set(_data))
    if len(_data) != 1 and len(_data) <= 10:
        colormap_data = get_iter_colormap(_data, 'tab10')
        for _d in _data:
            colors_data[_d] = next(colormap_data)
    else:
        colors_data[_data[0]] = color1
    # Prepare xtick-labels and colors

    if verbose is True:
        print('\nPreparing labels for: ')
    for s in SF:
        if exclude is not None:
            if type(exclude) == str:
                exclude = [exclude]
            SKIP = False
            for _mode in s.cst.keys():
                if _mode in exclude:
                    SKIP = True
            if SKIP is True:
                continue
        if verbose is True:
            print(s)
        name = s.stats.name.split('_')[0]
        # name =
        mname = None
        if s.stats.model in model:
            mname = name
            try:
                name = s.stats.model   # s.stats.name.split('_')[1]
            except Exception:
                name = s.stats.name

        if mode == 'cc':
            labels = s.stats.modes_cc_in.names
        elif mode == 'stoneley':
            labels = []
            for _mm in s.stats.modes_in:
                if _mm.sens.lower() == mode.lower():
                    labels += [_mm.name]
        else:
            labels = s.stats.modes_in.names

        if verbose is True:
            msg = 'Reading labels: %s' % labels
            print(msg)

        m_list = labels[:]  # removing modes from branch
        for m in m_list:

            # Check if label actually present in cst, else remove
            if m not in s.cst.keys():
                labels.remove(m)
                continue

            if str(degree) not in s.cst[m].keys():
                labels.remove(m)
                continue

            mm = split_digit_nondigit(m)

            if mode == 'cc':
                if l != 'all' and isinstance(l, (list,)):
                    if int(mm[2]) is l[0] and int(mm[6]) is l[1]:
                        pass
                    elif int(mm[2]) is l[1] and int(mm[6]) is l[0]:
                        pass
                    else:
                        labels.remove(m)
            else:
                if l != 'all' and isinstance(l, (int,)):
                    if int(mm[2]) is not l:
                        labels.remove(m)
                if n != 'all':
                    if type(n) is int:
                        if int(m[0]) is not n:
                            labels.remove(m)
                    else:
                        if int(m[0]) not in n:
                            labels.remove(m)
        labels = [format_name(x, 2) for x in labels]

        # Loops through all labels for each splitting function s
        # the labels are created before, using the information in
        # s.stats.modes_in.names
        # If label is not in x_ticklabel it will be append, if it is
        # in the list the corresponding model is assigned a marker
        # and color
        for label in labels:
            # print("label:", label, x_ticklabel)
            # print("SF:", s)

            if label not in x_ticklabel:
                if mode in ['s', 't', 'S', 'T']:
                    if mode.upper() not in label:
                        continue
                x_ticklabel.append(label)
                s_nums[label] = i
                i += 1

            # Why was this if condition here?
            # else:
            mlabel = s.stats.model
            # print("Model:", s, mlabel)
            # Sometimes data sneaks into the model plot? have to check
            # if the s.stats.model is not None here, else continue
            if mlabel is None:
                continue
            if (
                 s.stats.model not in label_model_set and
                 s.stats.model is not None
                 ):
                label_model_set[s.stats.model] = False
            if verbose is True:
                print('Labels ', label, mlabel, mark, markers)
            if mlabel not in mark:
                mark[mlabel] = markers[mlabel]
            if mlabel not in colors:
                if cmap == 'grey':
                    colors[mlabel] = 'grey'
                else:
                    try:
                        colors[mlabel] = next(colormap)
                    except StopIteration:
                        if label2 is not None:
                            colormap = get_iter_colormap(model[0:-1], cmap,
                                                         random_values='center')
                        else:
                            colormap = get_iter_colormap(model, cmap,
                                                         random_values='center')
                        colors[mlabel] = next(colormap)
            if model_colors is not None:
                if mlabel in model_colors:
                    colors[mlabel] = model_colors[mlabel]
    # Ordering modes
    _out = sort_modes(mode, x_ticklabel, ordering, verbose)
    x_ticklabel, labels, s_nums, mode_list = _out[:]

    # Prepare Vlines that separates overtones here
    set_line = get_vlines(mode, ordering, s_nums, labels, mode_list, verbose)
    x_init = range(len(s_nums))
    label_data_set = {}
    # Plotting loop
    if verbose is True:
        print('\nPlotting:')

    plot_dict = {}
    for s in SF:
        if s.stats.name not in label_data_set:
            label_data_set[s.stats.name] = False
        if kind == "cst":
            # cst
            sf = s.cst
            sf_err = s.cst_errors
        elif kind == 'dst':
            # dst
            sf = s.dst
            sf_err = s.dst_errors

            if len(sf_err) > 0:
                tmp_e = sf_err[list(sf_err.keys())[0]]['0']['uncertainty']
                if tmp_e == 0:
                    sf_err = s.cst_errors
        # Loop over mode names of csts, key=modename
        for key in sf.keys():
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
                if l != 'all' and isinstance(l, (int,)):
                    if int(m[2]) is not l:
                        skip = True
                if n != 'all':
                    if type(n) is int:
                        if int(m[0]) is not n:
                            skip = True
                    else:
                        if int(m[0]) not in n:
                            skip = True
            if mode == 'cc':
                m = split_digit_nondigit(key)
                if l != 'all' and isinstance(l, (list,)):
                    if int(m[2]) is l[0] and int(m[6]) is l[1]:
                        pass
                    elif int(m[2]) is l[1] and int(m[6]) is l[0]:
                        pass
                    else:
                        skip = True
            if mode.lower() == 'stoneley':
                for _mm in s.stats.modes_in:
                    if _mm.sens.lower() == mode.lower():
                        skip = False

            if skip is True:
                continue
            if str(degree) not in sf[key] and degree not in sf[key]:
                continue

            if order == 'all':
                coeffs = sf[key][str(degree)]
                try:
                    errors_temp = sf_err[key][str(degree)]
                    # print(s.stats.name, key, degree, errors_temp)
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

            elif order == 0:
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

            if parts != 'all':
                coeffs = [sf[key][str(parts[1])][parts[2]]]
            for i, (cst, err, erru, errl) in enumerate(zip(coeffs, errors,
                                                           errors_up,
                                                           errors_lw)):
                # eval = err
                # Check if degree is 0, then plot fc/Q
                if (degree == 0 and (plot_f_center is True or plot_Q is True)):
                    cst_fQ = None
                    if verbose:
                        print(key, type(key))
                    _mode = read_modes(modenames=str(key))[0]
                    if plot_f_center is True:
                        # fc = f0 + (4pi)**-1/2 * Re(c00) in microHerz
                        c00 = s.cst[key]['0']
                        fc = _mode.freq * 1e3 + 1./np.sqrt(4. * np.pi) * c00
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
                                msg = 'c00=%.2f %.2f/ %.2f'
                                print(msg % (c00, errl, erru))
                            erru = 1./np.sqrt(4. * np.pi) * abs(erru)
                            errl = 1./np.sqrt(4. * np.pi) * abs(errl)
                            if verbose is True:
                                print('f=%.2f -%.2f/ +%.2f\n' % (fc,
                                                                 errl, erru))
                            # eval = np.array([errl, erru]).reshape((2, 1))
                        if plot_Q is True:
                            cst_fQ = {'f': cst}
                    if plot_Q is True:
                        # print(s.stats.name, i, key, cst, err, erru, errl)
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
                                msg = 'd00=%.2f %.2f/ %.2f'
                                print(msg % (cst, errl, erru))

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
                            # eval = np.array([erru, errl]).reshape((2, 1))
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
                            weight = "bold"
                        elif plot_Q is True:
                            title = 'Quality factor'
                            weight = "bold"
                        elif kind == 'cst':
                            title = '$\mathrm{Re}(c_{%s0})$' % (degree)
                            weight = None
                        elif kind == 'dst':
                            title = '$\mathrm{Re}(d_{%s0})$' % (degree)
                            weight = None
                        elif parts != 'all':
                            title = '$\mathrm{%s}(c_{%s%s})$' % (parts[0], parts[1], parts[2])
                            weight = None
                        if subtitle is not None:
                            title += '\n{}'.format(subtitle)
                    else:
                        t = int(np.ceil(i/2.))
                        if i_col == 0:
                            if kind == 'cst':
                                title = '$\mathrm{Re}(c_{%s%s})$' % (degree, t)
                                weight = None
                            elif kind == 'dst':
                                title = '$\mathrm{Re}(d_{%s%s})$' % (degree, t)
                                weight = None
                        else:
                            if kind == 'cst':
                                title = '$\mathrm{Im}(c_{%s%s})$' % (degree, t)
                                weight = None
                            elif kind == 'dst':
                                title = '$\mathrm{Im}(d_{%s%s})$' % (degree, t)
                                weight = None

                    if isinstance(yaxs_split, list):
                        gridspec_kw = {'height_ratios': [1, 0.3]}
                        fig, (axt2, axt) = plt.subplots(
                                                2, 1, sharex=True,
                                                gridspec_kw=gridspec_kw,
                                                facecolor='w')
                        fig.set_size_inches(fig_size)
                    else:
                        axt = plt.subplot(gs[i_row, i_col])

                    if fig_abc:
                        if isinstance(yaxs_split, list):
                            axt2.set_title('%s) %s' % (fig_abc, title),
                                           weight=weight)
                        else:
                            axt.set_title('%s) %s' % (fig_abc, title),
                                          weight=weight)
                    else:
                        if isinstance(yaxs_split, list):
                            axt2.set_title(title, weight="bold")
                        else:
                            axt.set_title(title, weight="bold")

                    if plot_f_center is True:
                        if isinstance(yaxs_split, list):
                            axt2.set_ylabel('$\delta f$ ($\mu$Hz)')
                        else:
                            axt.set_ylabel('$\delta f$ ($\mu$Hz)')
                    elif plot_Q is True:
                        if isinstance(yaxs_split, list):
                            axt2.set_ylabel('$\delta Q$')
                        else:
                            axt.set_ylabel('$\delta Q$')
                    else:
                        if i_col == 0:
                            if isinstance(yaxs_split, list):
                                axt2.set_ylabel('$f$ ($\mu$Hz)')
                            else:
                                axt.set_ylabel('$f$ ($\mu$Hz)')
                    # axt.grid()
                    [i.set_linewidth(border_width) for i in axt.spines.values()]
                    axt.hlines(0, -1, len(s_nums), color='lightgray')
                    axt.set_xlim(-0.5, len(s_nums)-0.5)
                    axt.set_xticks(x_init)

                    axt.tick_params(labelbottom=False, labeltop=False,
                                    bottom=True, top=True,
                                    left=True, right=True)

                    # from IPythone import embed; embed()
                    # from matplotlib.ticker import StrMethodFormatter
                    # FMT = StrMethodFormatter("{:>4}")
                    if format_yticks is not None:
                        from matplotlib.ticker import FormatStrFormatter
                        FMT = FormatStrFormatter(format_yticks)
                        axt.yaxis.set_major_formatter(FMT)
                        if isinstance(yaxs_split, list):
                            axt2.yaxis.set_major_formatter(FMT)

                    if add_colorbar == 'bottom':
                        label_row = rows - 3
                    else:
                        label_row = rows - 1

                    if i_row == label_row or order == 0:
                        axt.tick_params(labelbottom=True,
                                        bottom=True, top=True,
                                        left=True, right=True)
                        axt.set_xticklabels(x_ticklabel, fontsize=12, zorder=0)
                        for tick in axt.get_xticklabels():
                            tick.set_rotation(rot)
                    elif i_row == i_col:
                        if labeltop is True:
                            axt.tick_params(labelbottom=False, labeltop=True,
                                            bottom=True, top=True,
                                            left=True, right=True)
                            axt.set_xticklabels(x_ticklabel, fontsize=12)
                            for tick in axt.get_xticklabels():
                                tick.set_rotation(rot)

                    if ylim != 'auto':
                        axt.set_ylim(ylim)

                    #  split y axis
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
                        axt.plot((-d, +d), (1-d, 1+d), **kwargs)  # bottom/left
                        axt.plot((1-d, 1+d), (1-d, 1+d), **kwargs)  # bottom/right

                        kwargs.update(transform=axt2.transAxes, color='k',
                                      clip_on=False)
                        axt2.plot((1-d, 1+d), (-d, +d), **kwargs)  # top/right
                        axt2.plot((-d, +d), (-d, +d), **kwargs)  # top/left

                        # axt = [axt, axt2]
                    ax_dict[i] = axt
                    plot_dict[axt] = {}

                # set Eval here
                if erru is not None and errorbars == 'asymmetrical':
                    eval = np.array([abs(errl), abs(erru)]).reshape((2, 1))
                elif errorbars == 'symmetrical':
                    eval = np.array([abs(err), abs(err)]).reshape((2, 1))
                if err is None:
                    eval = None

                ax = ax_dict[i]
                name = format_name(key, 2)
                mname = s.stats.model  # s.stats.name.split('_')[0]

                if mname in model and mname is not None:
                    if verbose is True:
                        print('current model: %s | %s' % (mname, name))
                        print(linestyle)
                    try:
                        s_num = s_nums[name]
                    except KeyError:
                        continue

                    # zorder = 10
                    # if mname == label2:
                    #     zorder = 150
                    # else:
                    #     eval = err

                    _label = s.stats.model
                    _width = mname

                else:
                    # color1 = colors_data[s.stats.name]
                    if verbose is True:
                        print('current data prep: %s, %i' % (name, i))
                    if name not in s_nums:
                        continue
                    s_num = s_nums[name]

                    # Check for current n of mode, if new branch starts
                    # draw a vertical line
                    if s_num in set_line:
                        ax.axvline(x=s_num - 0.5, color='lightgray',
                                   linestyle='--')
                        if isinstance(yaxs_split, list):
                            axt2.axvline(x=s_num - 0.5, color='lightgray',
                                         linestyle='--')

                    _label = s.stats.name
                    _width = dkey
                    if _label.lower() != 'data':
                        if _label.lower().startswith('data'):
                            _width = _label.split()[0]
                            if verbose:
                                print(_label, _label.split()[0], width[_width])
                            idx = _label.lower().split()[0].split('data')[-1]
                            if idx == '':
                                idx = 0
                            else:
                                idx = int(idx)
                            if idx + 1 > len(_label_damping):
                                _label_damping.append([])
                            if _label not in _label_damping[idx]:
                                _label_damping[idx] += [_label]

                _ax = ax_dict[i]
                if _label not in plot_dict[_ax]:
                    plot_dict[_ax][_label] = [[], [], [], []]
                plot_dict[_ax][_label][0].append(s_num + width[_width])
                try:
                    # f and Q plot are stored as arrays?
                    plot_dict[_ax][_label][1].append(cst[0])
                except IndexError:
                    # cst is a float
                    plot_dict[_ax][_label][1].append(cst)

                # Right now only symmetrical errorbars
                # Have to fix it for asymmetrical.. again
                try:
                    # f and Q plot are stored as arrays?
                    plot_dict[_ax][_label][2].append(err[0])
                except IndexError:
                    # err is a float
                    plot_dict[_ax][_label][2].append(err)
                plot_dict[_ax][_label][3].append(s)

    # Set the colors for the different dampings
    if len(_label_damping) != 0:
        cmap_damping = list([0] * len(_label_damping))
        zorder_damping = list([0] * len(_label_damping))
        for _j, _ld in enumerate(_label_damping):
            y = [float(y.split()[1]) for y in _ld]
            # x.sort()
            x = _ld
            _ld = [x for _, x in sorted(zip(y, x))]
            # _label_damping = ["data %s" % str(y) for y in x]
            if type(cmap_all_damping) != list:
                cad = cmap_all_damping
            else:
                cad = cmap_all_damping[_j]
            _cmap = getattr(cm, cad)
            _N = len(_ld)
            # had something to do with the dmapings?
            if skip2 is True:
                cmap_tmp = _cmap(np.linspace(0, 1, _N + 2))
                cmap_tmp = cmap_tmp[2:]
            else:
                cmap_tmp = _cmap(np.linspace(0, 1, _N))

            _cmap_damping = {}
            _zorder_damping = {}
            for _i, _l in enumerate(_ld):
                _cmap_damping[_l] = cmap_tmp[_i]
                _zorder_damping[_l] = 99 - _i
            cmap_damping[_j] = _cmap_damping
            zorder_damping[_j] = _zorder_damping

    legend_set = False
    # Do the plotting loop here using plot_dict
    # from IPython import embed; embed()
    for ax, data in plot_dict.items():
        for models, items in data.items():
            x = items[0]
            cst = items[1]
            eval = items[2]
            s = items[3]
            _linestyle = 'None'
            _markersize = markersize
            _label = models
            _plotstyle = 'errorbar'

            if models == 'data':
                _marker = marker
                _color = color1
                _zorder = 101
                _label = data_label[0]
            elif models.startswith('data'):
                # Trying to make it possible to plot all dampings for a run
                key = _label.lower().split()[0]
                idx = key.split('data')[-1]

                if idx == '':
                    idx = 0
                else:
                    idx = int(idx)
                _marker = markers[datain[key]]
                _color = cmap_damping[idx][_label]
                _zorder = zorder_damping[idx][_label]
                if _label.split()[-1] == '0':
                    _zorder = 100
                _label = datain[key]
                # if damping_label is not None:
                #     _label = label.split()[1]
                # else:
                #     _label = None
            else:
                _zorder = 10
                _marker = mark[models]
                # For Splitting functions with more then one mode it can happen
                # that there are S20/S40 values doubled in this array
                if models == 'S20RTS':
                    # Sorting the x and y values for the line plot here
                    cst = [y for _, y in sorted(zip(x, cst))]
                    x = sorted(x)
                    _plotstyle = 'plot'
                    _marker = 'None'
                    _linestyle = '-'
                    _linewidth = linewidth_S20
                    _color = 'grey'
                elif models == 'S40RTS':
                    # Sorting the x and y values for the line plot here
                    cst = [y for _, y in sorted(zip(x, cst))]
                    x = sorted(x)
                    _plotstyle = 'plot'
                    _linestyle = ':'
                    _marker = 'None'
                    _linewidth = linewidth_S20
                    _color = 'grey'
                elif models == 'SP12RTS':
                    # Sorting the x and y values for the line plot here
                    cst = [y for _, y in sorted(zip(x, cst))]
                    x = sorted(x)
                    _plotstyle = 'plot'
                    _linestyle = 'dashdot'
                    _marker = 'None'
                    _linewidth = linewidth_S20
                    _color = 'grey'
                else:
                    _zorder = 20
                    _color = colors[models]

            if _plotstyle == 'errorbar':
                # if isinstance(yaxs_split, list):
                #     ax, ax2 = ax[:]

                ax.errorbar(x, cst, yerr=eval,
                            marker=_marker,
                            markersize=_markersize, color=_color,
                            linestyle=_linestyle, elinewidth=elinewidth,
                            capsize=3, zorder=_zorder, label=_label)
                if isinstance(yaxs_split, list):
                    axt2.errorbar(x, cst, yerr=eval,
                                marker=_marker,
                                markersize=_markersize, color=_color,
                                linestyle=_linestyle, elinewidth=elinewidth,
                                capsize=3, zorder=_zorder, label=_label)
            else:
                ax.plot(x, cst,
                        marker=_marker,
                        markersize=_markersize, color=_color,
                        linestyle=_linestyle, linewidth=_linewidth,
                        zorder=_zorder, label=_label)
                if isinstance(yaxs_split, list):
                    axt2.plot(x, cst,
                            marker=_marker,
                            markersize=_markersize, color=_color,
                            linestyle=_linestyle, linewidth=_linewidth,
                            zorder=_zorder, label=_label)
        # removing whiskers from legend
        ehandles, elabels = ax.get_legend_handles_labels()
        # ehandles = [h[0] for h in ehandles]
        for i, h in enumerate(ehandles):
            try:
                ehandles[i] = h[0]
            except TypeError:
                # If ax.plot is used, h won't be subscriptable
                pass

        if N != 1 and order == 'all':
            if add_colorbar is not False:
                x_pos = cols - 1.9
            else:
                x_pos = cols - 0.9
            y_pos = 1
            if legend_show and legend_set is False:
                ax.legend(ehandles, elabels, frameon=False,
                          handletextpad=0.25,
                          bbox_to_anchor=(x_pos, y_pos),
                          loc=2, borderaxespad=0.)
        else:
            if legend_show and legend_set is False:
                ax.legend(ehandles, elabels, frameon=False,
                          handletextpad=0.25,
                          bbox_to_anchor=(1, 0.5),
                          loc='center left', borderaxespad=0.)
        legend_set = True

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=fs)
    if len(ax_dict) != 0:
        for ax in ax_dict.values():
            mpl.rcParams.update({'font.size': fs})
            ax.tick_params(axis='both', which='major', labelsize=fs)
            [i.set_linewidth(border_width) for i in ax.spines.values()]
            ax.xaxis.set_tick_params(which='both', width=tick_width)
            ax.yaxis.set_tick_params(which='both', width=tick_width)
            ax.xaxis.label.set_size(fs)
            ax.yaxis.label.set_size(fs)
            ax.title.set_size(fs)
    else:
        mpl.rcParams.update({'font.size': fs})
        ax.tick_params(axis='both', which='major', labelsize=fs)
        [i.set_linewidth(border_width) for i in ax.spines.values()]
        ax.xaxis.set_tick_params(which='both', width=tick_width)
        ax.yaxis.set_tick_params(which='both', width=tick_width)
        ax.title.set_size(fs)

    if isinstance(yaxs_split, list):
        mpl.rcParams.update({'font.size': fs})
        axt2.tick_params(axis='both', which='major', labelsize=fs)
        [i.set_linewidth(border_width) for i in axt2.spines.values()]
        axt2.xaxis.set_tick_params(which='both', width=tick_width)
        axt2.yaxis.set_tick_params(which='both', width=tick_width)
        axt2.title.set_size(fs)

    if add_colorbar:
        _cmap = getattr(cm, cmap_all_damping)
        cax = plt.subplot(gs[:, -1])
        cax.set_aspect(cbar_aspect)
        # Need to check how to properly get the colorvalues
        # plt.colorbar(_im, cax=cax, orientation='horizontal',
        #              label=right_y_label)

        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=_cmap,
                                norm=norm,
                                orientation='vertical')

        if type(colorbarlabel) is not str:
            multicolor_ylabel(cax, colorbarlabel,
                              colorbar_multicolor, 'y',
                              ylabel_anchor=colorbarlabel_anchor)

    if savefig is not False:
        if degree == 0 or order == 0:
            t = 0
        elif order == 'all':
            t = 't'
        fname = '%s_%s%s%s' % (mode, kind[0], degree, t)
        # plt.tight_layout()
        if savefig is True and type(savefig) is bool:
            fig.savefig('cbranch_%s.png' % fname, orientation='landscape',
                        dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches,
                        transparent=True)
        else:
            fig.savefig(**savefigargs, orientation='landscape',
                        dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches,
                        transparent=True)

    return SF


def get_vlines(mode, ordering, s_nums, labels, mode_list, verbose=False):
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
    return set_line


def sort_modes(mode, x_ticklabel, ordering, verbose=False):
    # Ordering modes
    mode_list = Modes()
    labels = []
    s_nums = {}
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
    return x_ticklabel, labels, s_nums, mode_list
