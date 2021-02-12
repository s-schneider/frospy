from frospy.core.splittingfunc import SplittingFunc as SplitF

# To Run on eejit, basemap is not installed yet
from frospy.core.splittingfunc.plot import _plot_coeffs, _plot_map, Bin
from frospy.core.splittingfunc.splittingfunc import _calc_SH_matrix
from frospy.core.splittingfunc.plot import sens_kernel, sensC_kernel

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
try:
    import pyshtools as sh
except ImportError:
    pass

from obspy.core import AttribDict
from frospy.util.base import sort_human
import os

import numpy as np
from copy import deepcopy

import fnmatch
import traceback


class Set(object):
    """
    Container for frospy.core.splittingfunc.SplittingFunc objects
    """
    def __init__(self, splitf=None):
        self.splitf = []
        if isinstance(splitf, SplitF):
            splitf = [splitf]
        if splitf is not None:
            self.splitf.extend(splitf)

    def __str__(self, extended=False):
        out = str(len(self.splitf)) + ' SplittingFunc(s) in Set:\n'
        if len(self.splitf) <= 20 or extended is True:
            out = out + "\n".join([_i.__str__() for _i in self])
        else:
            out = out + "\n" + self.splitf[0].__str__() + "\n" + \
                '...\n(%i other splitf)\n...\n' % (len(self.splitf) - 2)\
                + self.splitf[-1].__str__() + '\n\n[Use "print(' + \
                'Set.__str__(extended=True))" to print all Traces]'
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def __nonzero__(self):
        """
        A Set is considered zero if has no SplitFs.
        """
        return bool(len(self.splitf))

    def __len__(self):
        """
        Return the number of SplitFs in the Set object.

        .. rubric:: Example

        >>> Set = Set([SplitF(), SplitF(), SplitF()])
        >>> len(Set)
        3
        """
        return len(self.splitf)

    def __iter__(self):
        """
        Return a robust iterator for Set.splitf.

        Doing this it is safe to remove splitf from Sets inside of
        for-loops using Set's
        :meth:`~frospy.core.Set.Set.remove`
        method. Actually this creates a new iterator every time a splitf is
        removed inside the for-loop.

        """
        return list(self.splitf).__iter__()

    def __add__(self, other):
        """
        Add two Set or a Set with a single SplitF.

        :type other: :class:`~frospy.core.Set.Set` or
            :class:`~frospy.core.splitf.SplitF`
        :param other: Set or SplitF object to add.
        :rtype: :class:`~frospy.core.Sets.Set`
        :returns: New Set object containing references to the SplitFs of the
            original Set

        .. rubric:: Examples

        1. Adding two Set

            >>> st1 = Set([SplitF(), SplitF(), SplitF()])
            >>> len(st1)
            3
            >>> st2 = Set([SplitF(), SplitF()])
            >>> len(st2)
            2
            >>> Set = st1 + st2
            >>> len(Set)
            5

        2. Adding Set and SplitF

            >>> Set2 = st1 + SplitF()
            >>> len(Set2)
            4
        """
        if isinstance(other, SplitF):
            other = Set([other])
        if not isinstance(other, Set):
            raise TypeError
        splitf = self.splitf + other.splitf
        return self.__class__(splitf=splitf)

    def __iadd__(self, other):
        """
        Add two Sets with self += other.

        It will extend the current Stream object with the traces of the given
        Stream. Traces will not be copied but references to the original traces
        will be appended.

        :type other: :class:`~frospy.core.Set.Set` or
            :class:`~obspy.core.splitf.SplitF`
        :param other: Set or SplitF object to add.

        .. rubric:: Example

        >>> Set = Set([SplitF(), SplitF(), SplitF()])
        >>> len(Set)
        3

        >>> Set += Set([SplitF(), SplitF()])
        >>> len(Set)
        5

        >>> Set += SplitF()
        >>> len(Set)
        6
        """
        if isinstance(other, SplitF):
            other = Set([other])
        if not isinstance(other, Set):
            raise TypeError
        self.extend(other.splitf)
        return self

    def __getitem__(self, index):
        """
        __getitem__ method of frospy.Set objects.

        :return: SplitF objects
        """
        return self.splitf.__getitem__(index)

    def append(self, splitf):
        """
        Append a single SplitF object to the current Set object.
        """
        if isinstance(splitf, SplitF):
            self.splitf.append(splitf)
        else:
            msg = 'Append only supports a single SplittingFunc object as an '
            msg += 'argument.'
            raise TypeError(msg)
        return self

    def copy(self):
        return deepcopy(self)

    def del_duplicates(self):
        """
        Delete duplicate splitf in current Set object.
        Duplicate splitf have the same time window, frequency window,
        weight and station.
        """
        cop = self.copy()
        for splitf1 in self.__iter__():
            k = 0
            rem = Set()
            for splitf2 in cop.__iter__():
                if splitf1 == splitf2:
                    k = k + 1
                    rem.append(splitf2)
            if k >= 2:
                for b in range(1, k):
                    cop.remove(rem[b])
        return cop

    def extend(self, splitf_list):
        """
        Extend the current Set object with a list of SplitF objects.
        """
        if isinstance(splitf_list, list):
            for _i in splitf_list:
                # Make sure each item in the list is a trace.
                if not isinstance(_i, SplitF):
                    msg = 'Extend only accepts a list of SplitF objects.'
                    raise TypeError(msg)
            self.splitf.extend(splitf_list)
        elif isinstance(splitf_list, Set):
            self.splitf.extend(splitf_list.splitf)
        else:
            msg = 'Extend only supports a list of SplitF objects as argument.'
            raise TypeError(msg)
        return self

    def remove(self, splitf):
        """
        Remove the first occurrence of the specified SplitF object in the
        Set object. Passes on the remove() call to self.splitf.
        """
        self.splitf.remove(splitf)
        return self

    def select(self, damp=None, name=None, ncoeffs=None, nsegments=None,
               modes=None, modes_cc=None):
        """
        Return new Set object only with these splittingfunctions
        that match the given stats criteria
        (e.g. all splitf with ``damp="0.01"``).
        """
        SFs = []
        for sf in self:
            # skip trace if any given criterion is not matched
            if damp is not None:
                if not fnmatch.fnmatch(str(float(sf.stats.damp)),
                                       str(float(damp))):
                    continue
            if name is not None:
                if not fnmatch.fnmatch(sf.stats.name.upper(),
                                       name.upper()):
                    continue
            if ncoeffs is not None:
                if not fnmatch.fnmatch(str(float(sf.stats.ncoeffs)),
                                       str(float(ncoeffs))):
                    continue
            if nsegments is not None:
                if not fnmatch.fnmatch(str(float(sf.stats.nsegments)),
                                       str(float(nsegments))):
                    continue
            if modes is not None:
                i = 0
                for mname in sf.stats.modes_in.names:
                    if not fnmatch.fnmatch(mname.upper(), modes.upper()):
                        continue
                    else:
                        i += 1
                if i == 0:
                    continue
            if modes_cc is not None:
                i = 0
                for mname in sf.stats.modes_cc_in.names:
                    if not fnmatch.fnmatch(mname.upper(), modes.upper()):
                        continue
                    else:
                        i += 1
                if i == 0:
                    continue

            SFs.append(sf)
        return self.__class__(SFs)

    def set_mf(self, name, init_mf, final_mf):
        return

    def sort(self, keys=['damp', 'name', 'ncoeffs', 'nsegments', 'modes_in',
                         'modes_cc_in'], reverse=False):
        """
        Sort the picks in the Segment object.

        """
        # check if list
        msg = "keys must be a list of strings. Always available items to " + \
            "sort after: \n'damp', 'name', 'ncoeffs', 'nsegments', " + \
            "'modes_in', 'modes_cc_in'"
        if not isinstance(keys, list):
            raise TypeError(msg)

        # Loop over all keys in reversed order.
        for _i in keys[::-1]:
            self.splitf.sort(key=lambda x: getattr(x.stats, _i),
                             reverse=reverse)
        return self

    def plot(self, smin=None, smax=None, save_fig=False,
             order_by='degree', modelonly=True, cc_order='auto',
             **kwargs):
        """
        Plots coefficients of each self.splitf instance.
        Labels are set corresponding to self.splitf.stats.name

        order_by: 'degree' or 'coefficient'
        """

        if 'R' in kwargs:
            R = kwargs['R']
        else:
            R = -0.2

        if 'cmap' in kwargs:
            if hasattr(cm, kwargs['cmap']):
                cmap = getattr(cm, kwargs['cmap'])
            else:
                cmap = kwargs['cmap']
        else:
            if len(self) <= 9:
                cmap = getattr(cm, 'Set1')
            else:
                cmap = getattr(cm, 'tab20')

        if 'mode' in kwargs:
            plot_mode = kwargs['mode'].lower()
        else:
            plot_mode = 'all'
        # prepare colormap
        mcount = 0
        models = ['S20RTS', 'RR', 'REM', 'HT', 'AD', 'QM1', 'TZ']
        for splitf in self.splitf:
            if splitf.stats.name in models:
                mcount += 1
        # Create iterable colormap, will be used for data only
        # If model is input, grey dashed lines are used
        if type(cmap) is ListedColormap:
            if len(self) <= 9:
                cmap = cmap(np.linspace(0, 1, 9))
            else:
                cmap = cmap(np.linspace(0, 1, 20))

        elif type(cmap) is not str:
            if mcount != 0:
                cmap = cmap(np.linspace(0, 1, len(self)-mcount))
            else:
                cmap = cmap(np.linspace(0, 1, len(self)))
        # cmap_cst = iter(colormap)
        # cmap_dst = iter(colormap)

        # Initiate modes dict, contains all axes and figures
        modes_cst = AttribDict()
        modes_dst = AttribDict()

        # Set colormap counter
        cmap_count = 0
        for splitf in self.splitf:
            for mode, coeffs in splitf.cst.items():
                if plot_mode != 'all':
                    if mode.lower() != plot_mode.lower():
                        continue
                    if modelonly is False and cmap_count == 0:
                        if splitf.stats.name in models:
                            continue

                errors = splitf.cst_errors[mode]
                label = splitf.stats.name
                sdeg = sort_human(list(coeffs.keys()))[0]
                d00 = None
                d00_err = None

                # adjusting given smax to mode
                sdegs = [int(x) for x in coeffs.keys()]

                smax_in = smax
                if smax is None:
                    smax_in = max(sdegs)
                elif smax > max(sdegs):
                    smax_in = max(sdegs)

                if sdeg == '0':
                    if hasattr(splitf.dst, mode):
                        d00 = splitf.dst[mode]['0']
                        d00_err = splitf.dst_errors[mode]['0']

                if '-' in mode:
                    if mode.find('T') < mode.find('S') and cc_order == 'ST':
                        mode, coeffs = _swap_cc_order(mode, coeffs)
                    elif mode.find('T') > mode.find('S') and cc_order == 'TS':
                        mode, coeffs = _swap_cc_order(mode, coeffs)

                modes_cst = _plot_coeffs(
                            coeffs,  errors, mode, label, modes_cst, 'cst',
                            smin=smin, smax=smax_in,
                            colormap=cmap, colormap_index=cmap_count,
                            d00=d00, d00_err=d00_err, **kwargs)

            for mode, coeffs in splitf.dst.items():
                if plot_mode != 'all':
                    if mode.lower() != plot_mode.lower():
                        continue

                errors = splitf.dst_errors[mode]
                if splitf.stats.name == 'S20RTS':
                    label = "%s, R=%s" % (splitf.stats.name, R)
                else:
                    label = splitf.stats.name

                # d00 is already plotted in cst
                plot_coeffs = coeffs.copy()
                plot_coeffs.pop('0', None)
                plot_errors = errors.copy()
                plot_errors.pop('0', None)
                if len(plot_coeffs) == 0:
                    continue

                smin_in = smin
                if smin < 2 or smin is None:
                    smin_in = 2

                modes_dst = _plot_coeffs(
                            plot_coeffs, plot_errors, mode, label, modes_dst,
                            'dst', smin=smin_in, smax=smax,
                            colormap=cmap, colormap_index=cmap_count,
                            **kwargs)

            if splitf.stats.name not in models:
                cmap_count += 1

        if save_fig:
            for fi in plt.get_fignums():
                plt.figure(fi)
                filename = str(fi) + '_cst.png'
                fig = plt.gcf()
                fig.set_size_inches(12, 8)
                plt.savefig(filename, dpi=300, orientation='landscape')
                print("Saving figure %s" % filename)
            plt.close('all')

        return

    def plot_map(self, vmin=None, vmax=None, kind='cst',  # vlim='auto',
                 verbose=False, **kwargs):
        bins = Bin()
        smin = None
        smax = None

        if 'smin' in kwargs:
            if kwargs['smin'] != 'all':
                smin = int(kwargs['smin'])
        if 'smax' in kwargs:
            if kwargs['smax'] != 'all':
                smax = int(kwargs['smax'])

        if 'R' in kwargs:
            R = kwargs['R']

        if 'modes' in kwargs:
            modes = kwargs['modes']
            del kwargs['modes']

            if modes == 'cc':
                modes = []
                for sf in self:
                    for key in sf.cst.keys() and 'cst' in kind:
                        if '-' in key:
                            modes.append(key)
                    for key in sf.dst.keys() and 'dst' in kind:
                        if '-' in key:
                            modes.append(key)
            elif type(modes) != list:
                modes = [modes]

        elif kind == 'cst':
            modes = list(self.splitf[0].cst.keys())

        else:
            modes = list(self.splitf[0].dst.keys())

        if 'fig_abc' in kwargs:
            fig_abc = str(kwargs['fig_abc'])
            abc = ['-c', '-b', '-a', 'a', 'b', 'c',
                   'd', 'e', 'f', 'g', 'h', 'i']
            abc = iter(abc[abc.index(fig_abc)::])
            fig_abc = True
        else:
            fig_abc = False
        for i, mode in enumerate(modes):
            width = [0.8]  # kernel
            width.extend(np.ones(len(self.splitf))*3.2)
            width.append(0.15)  # colorbar
            gridspec_kw = {'width_ratios': width,
                           "left": 0.1, "bottom": 0.25, "right": 0.95,
                           "top": 0.95, "wspace": 0.05, "hspace": 0.}

            if vmin is None and vmax is None:  # and vlim == 'auto':
                vlim = [0, 0]
                for sf in self.splitf:
                    if (
                         kind == "cst" or
                         kind == "dst" and sf.stats.name != "AD"
                       ):
                        if verbose:
                            print('checking %s' % mode)
                        try:
                            if "-" in mode:
                                if len(sf.stats.modes_cc_in) > 1:
                                    m = sf.stats.modes_cc_in
                                    m = m.select(name=mode)[0]
                                else:
                                    m = sf.stats.modes_cc_in[0]
                            else:
                                m = sf.stats.modes_in.select(name=mode)[0]

                            if m.l != 0:
                                clm, sdegs = get_clm(sf, mode, kind,
                                                     smin, smax)
                                g, vminp, vmaxp = _convert_clm2g(clm)
                                if vminp < vlim[0]:
                                    vlim[0] = vminp
                                if vmaxp > vlim[1]:
                                    vlim[1] = vmaxp
                        except KeyError:
                            continue
                        if verbose:
                            print("vmin, vmax", vlim)
            else:
                vlim = [vmin, vmax]

            if os.path.exists(bins.sc_cstkernels):
                fig, axes = plt.subplots(nrows=1, ncols=len(self.splitf) + 2,
                                         gridspec_kw=gridspec_kw)
            else:
                fig, axes = plt.subplots(nrows=1, ncols=len(self.splitf) + 1)
                axes = np.hstack([np.array([None]), axes])

            for ax, sf in zip(axes.flat[1:-1], self.splitf):
                if (
                    kind == "cst" or
                    kind == "dst" and sf.stats.name != "AD"
                   ):
                    if verbose:
                        print('plotting %s' % mode)

                    if kind == 'dst' and sf.stats.name == 'S20RTS':
                        title = "%s, R=%s" % (sf.stats.name, R)
                    else:
                        title = sf.stats.name

                    if "-" in mode:
                        if len(sf.stats.modes_cc_in) > 1:
                            m = sf.stats.modes_cc_in.select(name=mode)[0]
                        else:
                            m = sf.stats.modes_cc_in[0]
                    else:
                        m = sf.stats.modes_in.select(name=mode)[0]

                    try:
                        if m.l != 0:
                            clm, sdegs = get_clm(sf, mode, kind, smin, smax)
                            if smax:
                                if smax > max(sdegs):
                                    kwargs['smax'] = max(sdegs)
                            else:
                                kwargs['smax'] = max(sdegs)

                            if smin:
                                if smin < 2:
                                    kwargs['smin'] = 2
                            else:
                                kwargs['smin'] = 2

                            im, fig = _plot_map(clm, m, kind=kind,
                                                suptitle=title,
                                                ax=ax, fig=fig,
                                                vmin=vlim[0], vmax=vlim[1],
                                                show_colorbar=False,
                                                show=False,
                                                **kwargs)
                            ax_cb = ax  # last plotted map axis saved
                        else:
                            ax.set_axis_off()
                    except Exception as e:
                        if verbose is True:
                            print(e)
                        ax.set_axis_off()
                else:
                    ax.set_axis_off()

            # sens_kernel plot

            ax = axes.flat[0]
            try:
                if m.l == 0 or max(sdegs) == 0 :
                    ax.set_axis_off()
                    if fig_abc and m.l != 0:
                        ax.set_title('%s)' % next(abc),
                                     x=0.2, y=1, weight="bold")
                elif "-" in mode:
                    if 'legend_show' not in kwargs:
                        kwargs['legend_show'] = True

                    if 'ticks' not in kwargs:
                        kwargs['ticks'] = False

                    try:
                        sensC_kernel(m, title=False, ax=ax, **kwargs)
                    except Exception:
                        e = traceback.format_exc()
                        print(e)

                    if kwargs['legend_show'] and fig_abc:
                        if 'fig_abc_y' in kwargs:
                            _y = kwargs['fig_abc_y']
                        else:
                            _y = 1.15
                        ax.set_title('%s)' % next(abc),
                                     x=0.2, y=_y, weight="bold")
                    elif fig_abc:
                        ax.set_title('%s)' % next(abc),
                                     x=0.2, y=1, weight="bold")

                # if "-" in mode or m.l == 0:
                #     ax.set_axis_off()
                #     if fig_abc and m.l != 0:
                #         ax.set_title('%s)' % next(abc),
                #                      x=0.2, y=1, weight="bold")
                else:
                    if 'legend_show' not in kwargs:
                        kwargs['legend_show'] = True

                    if 'ticks' not in kwargs:
                        kwargs['ticks'] = False
                    try:
                        sens_kernel(m, title=False, ax=ax, **kwargs)
                    except Exception:
                        e = traceback.format_exc()
                        print(e)

                    # I have to discuss with Su about this title
                    # I cannot see the title with y=1.15
                    if kwargs['legend_show'] and fig_abc:
                        if kwargs['fig_abc_y']:
                            _y = kwargs['fig_abc_y']
                        else:
                            _y = 1.15
                        ax.set_title('%s)' % next(abc),
                                     x=0.2, y=_y, weight="bold")
                    elif fig_abc:
                        ax.set_title('%s)' % next(abc),
                                     x=0.2, y=1, weight="bold")
                # colorbar position
                if m.l != 0:
                    # ax = ax_cb  # plot in the last map of each mode
                    # divider = make_axes_locatable(ax)
                    # cax = divider.append_axes("right", size="5%", pad=0.05)
                    cax = axes.flat[-1]
                    cax.set_aspect(12.5, adjustable='datalim')  # shrinking cb
                    if vlim[0] < 0 and vlim[1] > 0:
                        ticks = [vlim[0], 0., vlim[1]]
                    else:
                        ticks = [vlim[0], vlim[1]]
                    cb = fig.colorbar(im, ticks=ticks, format='%3.1f',
                                      ax=axes.ravel().tolist(), cax=cax)
                    cb.ax.set_title(r'$\mu$Hz', y=0.97)
                if m.l == 0:
                    cax = axes.flat[-1]
                    cax.set_axis_off()

                # figure size
                if 'xlen' in kwargs:
                    xlen = kwargs['xlen']
                elif len(self.splitf) == 1:
                    xlen = 3#5
                elif len(self.splitf) == 2:
                    xlen = 5#8
                elif len(self.splitf) == 3:
                    xlen = 7.5#12
                else:
                    xlen = 10#16

                if 'ylen' in kwargs:
                    ylen = kwargs['ylen']
                else:
                    ylen = 1.5

                if 'fig_size' in kwargs:
                    xlen, ylen = kwargs['fig_size']

                fig.set_size_inches(xlen, ylen)

            except ValueError as e:
                raise Exception(e)
            if 'savefigure' in kwargs and m.l != 0:
                save = kwargs['savefigure']
                if 'fileformat' in kwargs:
                    suffix = kwargs['fileformat']
                else:
                    suffix = 'png'

                if 'filename' in kwargs and save:
                    name = kwargs['filename']
                    fname = 'Set_%s_%s_%s' % (mode, kind, name)
                    fig.savefig('%s.%s' % (fname, suffix), bbox_inches="tight",
                                orientation='landscape', dpi=400,
                                pad_inches=0.01)
                elif save:
                    fname = 'Set_%s_%s' % (mode, kind)
                    fig.savefig('%s.%s' % (fname, suffix), bbox_inches="tight",
                                orientation='landscape', dpi=400,
                                pad_inches=0.01,)

        return fig


def get_clm(sf, mode, kind, smin, smax):
    if kind == 'cst':
        coeffs = sf.cst[mode]
    else:
        coeffs = sf.dst[mode]

    sdegs = [int(x) for x in coeffs.keys()]
    if kind == 'cst':
        ssum = sum([sum(x) for x in coeffs.values()])
        if max(sdegs) == 0 or ssum == 0.:
            raise Exception('smax=0')
    else:
        # skip degree 0
        if len(coeffs.keys()) in [0, 1]:
            raise Exception('smax=0')
        # dst Exception for cst inversion above 20
        for sdeg, cval in coeffs.iteritems():
            if int(sdeg) > 20:
                coeffs[sdeg] = np.zeros(2*int(sdeg)+1)

    SHmat = _calc_SH_matrix(coeffs, smin, smax)
    clm = sh.SHCoeffs.from_array(SHmat,
                                 normalization='ortho',
                                 csphase=-1)
    return clm, sdegs


def _convert_clm2g(clm):
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
    return g, vminp, vmaxp


def _swap_cc_order(mode, coeffs):
    mode = '-'.join(mode.split('-')[::-1])
    for d, c in coeffs.items():
        coeffs[d] = -c
    return mode, coeffs
