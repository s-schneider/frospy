# -*- coding: utf-8 -*-
"""
Module for handling nmPy SplittingFunc objects.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import, print_function
from copy import deepcopy

from obspy.core.util import AttribDict
try:
    import pyshtools as sh
except ImportError:
    pass

from frospy.core.splittingfunc.plot import _plot_coeffs, _plot_map
import matplotlib.pyplot as plt
from frospy.core.segment import Segment
from frospy.core.modes import Modes, Mode, format_name
from frospy.util.read import (read_modes_in, read_lines)
from frospy.util.write import write_pickle
from frospy.util.nmfile import sort_AttribDict
from frospy.util.base import (
            sc_degrees, cc_degrees, max_sc_degrees, sort_human,
            N_splitting_coeffs
            )

from frospy.core.splittingfunc.read import (read_cst, get_cst)
from obspy.core.event.base import QuantityError
import re
import numpy as np
import os
from collections import OrderedDict
import sys
import glob
import json


class Stats(AttribDict):
    defaults = {
        'name': None,
        'damp': None,
        'nsegments': 0,
        'ncoeffs': 0,
        'model': None,
        'modes_in': Modes(),
        'modes_out': Modes(),
        'modes_cc_in': Modes(),
        'modes_cc_out': Modes(),
        'initial_mf': None,
        'final_mf': None
    }

    _refresh_keys = {'modes_in', 'modes_out', 'modes_cc_in', 'modes_cc_out',
                     'nsegments', 'ncoeffs', 'damp', 'name', 'model'}

    def __init__(self, header={}):
        """
        """
        # super() is used to initialize AttribDict within Stats(),
        # which sets the defaults as initial dictionary
        # equivalent to AttribDict.__init__(self, header)
        super(Stats, self).__init__(header)

    def __setitem__(self, key, value):
        """
        """
        if key in self._refresh_keys:
            if key in ['damp', 'ncoeffs', 'nsegments']:
                value = value
            elif key in ['modes_in', 'modes_out', 'modes_cc_in',
                         'modes_cc_out']:
                if isinstance(value, Modes) or isinstance(value, Mode):
                    value = value
                elif type(value) is list:
                    value = [str(val) for val in value]
                else:
                    value = str(value)
            elif key == 'segment':
                if isinstance(value, Segment):
                    value = value
                elif type(value) is list:
                    value = [str(val) for val in value]
                else:
                    value = str(value)
            # equivalent to AttribDict.__setitem__(self, key, value)
            super(Stats, self).__setitem__(key, value)

    __setattr__ = __setitem__

    def __str__(self, extended=False):
        """
        Return better readable string representation of Stats object
        """
        _pretty_str = 'Damping: %s | Seg-Count: %s | ncoeff: %s \n'
        _pretty_str = _pretty_str % (self.damp, self.nsegments,
                                     self.ncoeffs)

        if isinstance(self.modes_in, Modes):
            _pretty_str += 'Mode(s):'
            for mod in self.modes_in:
                _pretty_str += ' %s' % (mod.name)

        if isinstance(self.modes_cc_in, Modes):
            _pretty_str += '\nCC mode(s):'
            for mod in self.modes_cc_in:
                _pretty_str += ' %s' % (mod.name)
        _pretty_str += '\nname: %s\n' % self.name

        return _pretty_str

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))


class SplittingFunc(object):
    defaults = {
        'cst': AttribDict(),
        'dst': AttribDict(),
        'cst_errors': AttribDict(),
        'dst_errors': AttribDict()
    }

    _refresh_keys = {'cst', 'dst', 'cst_errors', 'dst_errors'}

    def __init__(self, header=None, cst=None, dst=None,
                 cst_errors=None, dst_errors=None):
        if cst is None:
            cst = {}
        if dst is None:
            dst = {}

        self.cst = AttribDict(cst)
        self.dst = AttribDict(dst)

        if cst_errors is None:
            # If not given, set errors to 0, for plotting
            cst_errors = self.cst.copy()
            for m, degs in cst_errors.items():
                for deg, c in degs.items():
                    err = QuantityError()
                    err.uncertainty = np.zeros(len(c))
                    err.upper_uncertainty = np.zeros(len(c))
                    err.lower_uncertainty = np.zeros(len(c))
                    err.confidence_level = 0
                    cst_errors[m][deg] = err
        if dst_errors is None:
            dst_errors = self.dst.copy()
            for m, degs in dst_errors.items():
                for deg, c in degs.items():
                    err = QuantityError()
                    err.uncertainty = np.zeros(len(c))
                    err.upper_uncertainty = np.zeros(len(c))
                    err.lower_uncertainty = np.zeros(len(c))
                    err.confidence_level = 0
                    dst_errors[m][deg] = err

        self.cst_errors = AttribDict(cst_errors)
        self.dst_errors = AttribDict(dst_errors)

        if header is None:
            header = {}
        header = deepcopy(header)
        header.setdefault('ncoeffs', len_coeff(cst, dst))
        self.stats = Stats(header)

    def __str__(self, extended=False):
        _pstr = ''
        if self.stats.name is not None:
            _pstr += "%s with " % self.stats.name

        _pstr += "cst/dst for: "

        if isinstance(self.stats.modes_in, Modes):
            for mod in self.stats.modes_in:
                if mod.name in self.cst or mod.name in self.dst:
                    _pstr += " %s" % mod.name

        if isinstance(self.stats.modes_cc_in, Modes):
            for mod in self.stats.modes_cc_in:
                if mod.name in self.cst or mod.name in self.dst:
                    _pstr += " %s" % mod.name

        return _pstr

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def copy(self):
        return deepcopy(self)

    def write(self, filename, overwrite=False, format='pickle',
              file_per_mode=False):
        """
        writes cst file, order of modes is the same as order in
        self.cst.keys()
        """
        def write_dat(fname):
            with open(fname, 'w') as fh:
                for val in cf:
                    fh.write("{}\n".format(val))
            return

        try:
            # Get filename
            if not overwrite:
                i = 0
                while True:
                    msg = ''
                    path = os.path.dirname(filename)
                    f = os.path.basename(filename)
                    f, ext = os.path.splitext(f)

                    if os.path.exists(filename):
                        msg += '\033[93mSplittingfunction-file exist,'
                        msg += 'not overwriting\033[0m'
                        if i == 0:
                            f = "%s_%s" % (f, str(i))
                        i += 1
                        a = "_%s" % str(i-1)
                        b = "_%s" % str(i)
                        f = f.replace(a, b)
                        filename = os.path.join(path, "%s%s" % (f, ext))
                    else:
                        print(msg)
                        break

            if format == 'pickle':
                write_pickle(self, filename)

            elif format == 'json':
                data = {'SF': {}}
                for m, cst in self.cst.items():
                    data['SF'][m] = {}
                    for deg, c in cst.items():
                        data['SF'][m][deg] = c.tolist()
                with open(filename, 'w') as fh:
                    json.dump(data, fh)

            elif format == 'dat':
                if file_per_mode is False:
                    cf = []
                    for mode, c in self.cst.items():
                        degs = sorted(list(self.cst[mode].keys()))
                        for deg in degs:
                            coeff = self.cst[mode][deg]
                            if deg == '0':
                                cf += [coeff[0]]
                                cf += [self.dst[mode]['0'][0]]
                            else:
                                cf += coeff.tolist()
                    write_dat(filename)

                else:
                    for mode, c in self.cst.items():
                        cf = []
                        degs = sorted(list(self.cst[mode].keys()))
                        for deg in degs:
                            coeff = self.cst[mode][deg]
                            if deg == '0':
                                cf += [coeff[0]]
                                cf += [self.dst[mode]['0'][0]]
                            else:
                                cf += coeff.tolist()
                        print(mode)
                        write_dat("{}.dat".format(mode))

            else:
                raise IOError('Only support for pickle and json files.')

        except IOError:
            msg = "\033[91mCan't save file\n"
            msg += "Error message: %s\033[0m" % sys.exc_info()[1]
            print(msg)
        return

    def plot_map(self, fs=10, cmap='lies', **kwargs):
        smin = None
        smax = None
        show = True
        figs = []
        try:
            title = '%s, d=%s' % (int(self.stats.name), self.stats.damp)
        except Exception:
            title = self.stats.name

        if 'smin' in kwargs:
            if kwargs['smin'] != 'all':
                smin = int(kwargs['smin'])
        if 'smax' in kwargs:
            if kwargs['smax'] != 'all':
                smax = int(kwargs['smax'])
        if 'show' in kwargs:
            show = kwargs['show']
        if 'modes' in kwargs:
            modes = kwargs['modes']
            if modes == 'cc':
                modes_cst = []
                modes_dst = []
                for key in self.cst.keys():
                    if '-' in key:
                        modes_cst.append(key)
                for key in self.dst.keys():
                    if '-' in key:
                        modes_dst.append(key)

            elif type(modes) is not list:
                modes = modes
                modes_cst = modes_dst = modes

            del kwargs['modes']

        else:
            modes_cst = self.cst.keys()
            modes_dst = self.dst.keys()

        try:  # cst
            im_cst = []
            for mode in modes_cst:
                coeffs = self.cst[mode]
                if "-" in mode:
                    if len(self.stats.modes_cc_in) > 1:
                        m = self.stats.modes_cc_in.select(name=mode)[0]
                    else:
                        m = self.stats.modes_cc_in[0]
                else:
                    m = self.stats.modes_in.select(name=mode)[0]

                sdegs = [int(x) for x in coeffs.keys()]
                ssum = sum([sum(x) for x in coeffs.values()])
                if max(sdegs) == 0 or ssum == 0.:
                    continue
                SHmat = _calc_SH_matrix(coeffs, smin, smax)

                clm = sh.SHCoeffs.from_array(SHmat, normalization='ortho',
                                             csphase=-1)
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

                im, fig = _plot_map(clm, m, kind='cst', suptitle=title, fs=fs,
                                    cmap=cmap, **kwargs)
                im_cst.append(im)
                figs.append(fig)

        except (TypeError, KeyError):
            im = None

        try:  # dst
            im_dst = []
            for mode in modes_dst:
                coeffs = self.dst[mode]
                if "-" in mode:
                    if len(self.stats.modes_cc_in) > 1:
                        m = self.stats.modes_cc_in.select(name=mode)[0]
                    else:
                        m = self.stats.modes_cc_in[0]

                sdegs = [int(x) for x in coeffs.keys()]
                # skip degree 0
                if len(coeffs.keys()) in [0, 1]:
                    continue
                SHmat = _calc_SH_matrix(coeffs, smin, smax)

                clm = sh.SHCoeffs.from_array(SHmat, normalization='ortho',
                                             csphase=-1)
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

                im, fig = _plot_map(clm, m, kind='dst', suptitle=title,
                                    **kwargs)
                im_dst.append(im)
                figs.append(fig)

        except (TypeError, KeyError):
            im = None

        if show:
            plt.show()
        return im_cst, im_dst, figs

    def plot(self, savefig=False, smin=None, smax=None, colormap='red',
             return_ax=False,
             **kwargs):

        modes_cst = AttribDict()
        modes_dst = AttribDict()
        for mode, coeffs in self.cst.items():

            errors = self.cst_errors[mode]
            # There is a mode, coeffs mixup here, example 0T5-1S3
            label = self.stats.name
            sdeg = sort_human(list(coeffs.keys()))[0]
            d00 = None
            d00_err = None

            # adjusting given smax to mode
            sdegs = [int(x) for x in coeffs.keys()]

            smax_in = smax
            if smax is not None:
                if smax > max(sdegs):
                    smax_in = max(sdegs)
            else:
                smax_in = max(sdegs)

            if sdeg == '0':
                if hasattr(self.dst, mode):
                    if sdeg in self.dst[mode]:
                        d00 = self.dst[mode]['0']
                        d00_err = self.dst_errors[mode]['0']
            modes_cst = _plot_coeffs(
                        coeffs, errors, mode, label, modes_cst, 'cst',
                        d00=d00, d00_err=d00_err, smin=smin, smax=smax_in,
                        colormap=colormap, **kwargs
                        )

        for mode, coeffs in self.dst.items():

            errors = self.dst_errors[mode]
            # d00 is already plotted in cst
            label = self.stats.name
            plot_coeffs = coeffs.copy()
            plot_coeffs.pop('0', None)
            plot_errors = errors.copy()
            plot_errors.pop('0', None)
            if len(plot_coeffs) == 0:
                continue

            smin_in = smin
            if smin is None:
                if smin > 2:
                    smin_in = 2

            modes_dst = _plot_coeffs(
                        plot_coeffs, plot_errors, mode, label, modes_dst,
                        'dst', smin=smin_in, smax=smax, colormap=colormap,
                        **kwargs
                        )

        if savefig:
            for fi in plt.get_fignums():
                plt.figure(fi)
                fig = plt.gcf()
                fname = 'coeff_ReIm_%s' % (fi)
                fig.set_size_inches(12, 8)
                fig.savefig('%s.ps' % fname, orientation='landscape',
                            bbox_inches='tight')
        else:
            plt.show()

        if return_ax is True:
            return modes_cst, modes_dst
        else:
            return

    def get_fQ(self, mode_name):
        """
        Calculates center frequency and Q value of a mode in
        splitting function / self.

        return: (f_center, f err, f upper error, f lower error)
                (Q, Q error, Q upper error, Q lower error)
        """
        c00 = self.cst[mode_name]['0']
        err = self.cst_errors[mode_name]['0']['uncertainty']
        erru = self.cst_errors[mode_name]['0']['upper_uncertainty']
        errl = self.cst_errors[mode_name]['0']['lower_uncertainty']

        d00 = self.dst[mode_name]['0']
        derr = self.cst_errors[mode_name]['0']['uncertainty']
        derru = self.cst_errors[mode_name]['0']['upper_uncertainty']
        derrl = self.cst_errors[mode_name]['0']['lower_uncertainty']

        mode = self.stats.modes_in.select(name=mode_name)[0]
        fc = mode.freq * 1e3 + 1. / np.sqrt(4. * np.pi) * c00
        if err is not None:
            err = 1. / np.sqrt(4. * np.pi) * err

        if erru is not None:
            erru = 1. / np.sqrt(4. * np.pi) * abs(erru)
            errl = 1. / np.sqrt(4. * np.pi) * abs(errl)

        Qc = calc_Q(mode, fc, d00)

        if derr is not None:
            Qe = calc_Q(mode, fc, d00, derr)
            derr = Qc - Qe

        if derru is not None and derrl is not None:
            ## err in d00 translates to opposite sign error in Q
            #if derrl > 0:
            #    # if its positive only upper uncsrt exist
            #    derrl = 0
            #if derru < 0:
            #    # if its negative only lower uncsrt exist
            #    derru = 0

            # errl/erru they are flipped on purpose
            Qe_l = calc_Q(mode, fc, d00, derru)
            Qe_u = calc_Q(mode, fc, d00, derrl)
            derru = abs(Qc - Qe_u)
            derrl = abs(Qc - Qe_l)

        return (fc, err, erru, errl), (Qc, derr, derru, derrl)


def _calc_SH_matrix(coeffs, smin, smax):

    lsize = int(max(np.array(list(coeffs.keys()), dtype=int))) + 1
    SHmat = np.zeros([2, lsize, lsize])

    for sdeg, cval in coeffs.items():
        sdeg = int(sdeg)
        cval = iter(cval)
        # Exclude degree 0
        for mval in range(sdeg+1):
            if mval == 0 and sdeg == 0:
                SHmat[0, sdeg, mval] = 0.
#                _ = next(cval)
            elif mval == 0:
                SHmat[0, sdeg, mval] = next(cval)
            else:
                SHmat[0, sdeg, mval] = next(cval)
                SHmat[1, sdeg, mval] = -next(cval)

    if smin is not None:
        SHmat[0, 1:smin, :] = 0.
        SHmat[1, 1:smin, :] = 0.
    if smax is not None:
        SHmat[0, smax+1:, :] = 0.
        SHmat[1, smax+1:, :] = 0.
    return SHmat


def get_model_cst(setup, model, mode, smin, smax, kind="cst"):
    if model == 'PREM':
        if smin == 0:
            N = N_splitting_coeffs(smax, 0) + 1
        else:
            N = N_splitting_coeffs(smax, smin)
        return np.zeros(N)

    elif model in ('S20RTS', 'AD', 'RR'):
        mode = format_name(mode)
        format = model
        cst_out = read_cst(setup=setup, cfile=format, R=setup.dst_R)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        c = np.array([])
        if kind == "cst":
            for i in range(smin, smax+1, 2):
                if i == 0:
                    if model == 'RR':
                        c = np.hstack((0, 0))
                    else:
                        c = np.hstack((cst[mode][str(i)], dst[mode][str(i)]))
                else:
                    if model == 'RR' and str(i) not in cst[mode]:
                        c = np.hstack((c, np.zeros(2*i + 1)))
                    else:
                        c = np.hstack((c, cst[mode][str(i)]))
        if kind == "dst":
            for i in range(smin, smax+1, 2):
                c = np.hstack((c, dst[mode][str(i)]))
    elif type(model) == list:
        if smin == 0:
            N = N_splitting_coeffs(smax, 0) + 1
        else:
            N = N_splitting_coeffs(smax, smin)
        if len(model) == N:
            c = model
        else:
            raise IOError('expected model size: %s; received model size' % (N, len(c)))
    return c


def get_header(dir, modes_sc, modes_cc, name=None, damp=None, model=None):

    if dir is not None:
        # Find number of segments used
        events = glob.glob('%s/???????.misfit' % dir)
        seg_cnt = 0
        for file in events:
            seg_cnt += int(len(read_lines(file)) / 4.)
    else:
        events = None
        seg_cnt = None

    if modes_cc is None:
        modes_cc = Modes([])

    header = {'modes_in': modes_sc, 'modes_cc_in': modes_cc, 'damp': damp,
              'nsegments': seg_cnt, 'name': name, 'model': model}

    return header


def len_coeff(cst, dst):
    N = 0
    for mode in cst.values():
        for coeff in mode.values():
            if type(coeff) == np.ndarray:
                N += len(coeff)
            else:
                N += 1

    for mode in dst.values():
        for coeff in mode.values():
            if type(coeff) == np.ndarray:
                N += len(coeff)
            else:
                N += 1

    return N


def get_cst_list(infile, modes_dir, degree):
    cst_list = []
    dst_list = []
    for fh, md in zip(infile, modes_dir):
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = read_cst(fh, md)
        cst_sorted = AttribDict()
        dst_sorted = AttribDict()

        for k in cst.keys():
            cst_tmp = sort_AttribDict(cst[k], 'int')
            cst_sorted[k] = cst_tmp
        cst_list.append(cst_sorted)

        for k in dst.keys():
            dst_tmp = sort_AttribDict(dst[k], 'int')
            dst_sorted[k] = dst_tmp
        dst_list.append(dst_sorted)

    for cst in cst_list:
        for name, coeff in cst.items():
            new_coeff = []
            for deg in coeff:
                if degree == 'all':
                    new_coeff.append(deg)
                elif int(list(deg.keys())[0]) in degree:
                    new_coeff.append(deg)
            cst[name] = new_coeff
    for dst in dst_list:
        for name, coeff in dst.items():
            new_coeff = []
            for deg in coeff:
                if degree == 'all':
                    new_coeff.append(deg)
                elif int(list(deg.keys())[0]) in degree:
                    new_coeff.append(deg)
            dst[name] = new_coeff

    return cst, dst


def get_cst_d00_list(infile, modes_dir, degree):

    cst_list = []
    d00_list = []
    models = ['RR', 'TZ', 'S20RTS', 'REM']
    for fh, md in zip(infile, modes_dir):
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = read_cst(fh, md)
        cst_sorted = AttribDict()
        d00_sorted = AttribDict()
        for k in cst.keys():
            cst_tmp = sort_AttribDict(cst[k], 'int')
            cst_sorted[k] = cst_tmp
            if k.find("-") == -1 and fh not in models:
                d00_tmp = sort_AttribDict(dst[k], 'int')
                d00_sorted[k] = d00_tmp
        cst_list.append(cst_sorted)
        d00_list.append(d00_sorted)

    if degree != 'all':
        for cst in cst_list:
            for name, coeff in cst.items():
                new_coeff = []
                for deg in coeff:
                    if int(list(deg.keys())[0]) in degree:
                        new_coeff.append(deg)
                cst[name] = new_coeff

        for dst in d00_list:
            for name, coeff in dst.items():
                new_coeff = []
                for deg in coeff:
                    if int(list(deg.keys())[0]) in degree:
                        new_coeff.append(deg)
                dst[name] = new_coeff
    return cst_list, d00_list


def calc_cst_indices(modes_dir):
    """
    Returns a list of all end indices for a given modes.in and modes_cc.in file
    """

    def print_ind_names(inds, names, i_start, c_list, verbose=False):
        nami = iter(names)
        for i in inds:
            name = next(nami)

            if '-' in name:
                _name = name.split('-')
                mode_in = []
                for m in _name:
                    mode_in += re.split('(\D+)', m)[:]
                mode_in.append(i)
                mode_in.append(0)

                c = cc_degrees(mode_in)
                ind = 0
                for _c in c:
                    ind += 2 * _c + 1

            else:
                ind = sc_degrees(i)
                if ind == 0:
                    continue
                ind -= 1

            i_start += ind
            if verbose:
                print("%i, %s" % (i_start, name))
            c_list.append((i_start, name))

        return i_start, c_list

    m_sc, m_cc = read_modes_in(modes_dir)

    cst_inds = []
    cst_cc_inds = []
    dst_inds = []
    dst_cc_inds = []
    names = []
    names_cc = []
    for m in m_sc[1:]:
        mode = m.split()
        names.append(''.join(mode[:3]))

        cst_max = int(mode[3])
        mode_smax = max(max_sc_degrees(int(cst_max)))
        if cst_max > mode_smax:
            cst_max = mode_smax
        cst_inds.append(cst_max)

        dst_max = int(mode[4])
        if dst_max > mode_smax:
            dst_max = mode_smax
        dst_inds.append(dst_max)

    if m_cc is not None:
        no_cc = int(m_cc[0])
        for m in m_cc[1:no_cc + 1]:
            mode = m.split()
            names_cc.append(''.join(mode[0:3]) + '-' + ''.join(mode[3:6]))

            cst_cc_max = int(mode[6])
            cst_cc_inds.append(cst_cc_max)

            dst_cc_max = int(mode[7])
            dst_cc_inds.append(dst_cc_max)
    coef_order = []
    ind, coef_order = print_ind_names(cst_inds, names, 0, coef_order)
    if m_cc is not None:
        ind, coeff_order = print_ind_names(cst_cc_inds, names_cc, ind,
                                           coef_order)
    ind, coef_order = print_ind_names(dst_inds, names, ind, coef_order)

    if m_cc is not None:
        ind, coef_order = print_ind_names(dst_cc_inds, names_cc, ind,
                                          coef_order)

    return coef_order


def cst4cstmap(infile, modes_dir):
    if type(infile) != list:
        infile = [infile]
    cst = get_cst_list(infile, modes_dir, 'all')[0]

    for key in cst.keys():
        filename = "cst_%s.dat" % key
        with open(filename, 'w') as fh:
            if '-' in key:
                fh.write('0.0\n')
            for i, deg in enumerate(cst[key]):
                if i == 1:
                    fh.write('0.0\n')
                for d, csts in deg.items():
                    for value in csts:
                        line = "%e\n" % value
                        fh.write(line)
    return


def get_covariances(covmatrix, modes_dir):

    cst_cov_list = []
    if type(covmatrix) != list:
        covmatrix = [covmatrix]

    for covfile, md in zip(covmatrix, modes_dir):
        if covfile is None:
            cst_cov_list.append(None)
            continue
        cov = np.fromfile(covfile, dtype=np.dtype('f8'))
        M = int(np.sqrt(cov.shape[0]))
        cov = np.sqrt(cov.reshape(M, M).diagonal())
        # Read mode files
        modes, modes_cc = read_modes_in(md)
        cst_cov, dst_cov = get_cst(modes, modes_cc, cov, noc=True)

        cst_sorted = AttribDict()

        for k in cst_cov.keys():
            cst_tmp = sort_AttribDict(cst_cov[k], 'int')
            cst_sorted[k] = cst_tmp
        cst_cov_list.append(cst_sorted)

    return cst_cov_list, None


def build_deg_dependent_damp(setup, modes='all', degs='all', function='linear',
                             smag_diff=10):
    def get_m(x_min, x_max, function, x, y_max=smag_diff - 1):
        if x_max == x_min:
            return 0.
        if function == 'linear':
            a = y_max / x_max
            return 1. + x * a
        elif function == 'quad':
            a = y_max / x_max ** 2.
            return 1. + a * x**2.
        elif function == 'sqrt':
            a = y_max / np.sqrt(x_max)
            return 1. + a * np.sqrt(x)
        elif function == 'exp':
            a = y_max / np.log(x_max + 1)
            return 1. + a * np.log(x + 1)
        else:
            return function

    if modes != 'all':
        if type(modes) is not list:
            modes = [modes]
        modes = [m.lower() for m in modes]

    if degs != 'all':
        if type(degs) != list:
            degs = [int(degs)]
        else:
            degs = [int(x) for x in degs]

    mdamp = OrderedDict()
    for kind in setup.mzero.keys():
        mdamp[kind] = OrderedDict()
        for mode in setup.mzero[kind]:
            mdamp[kind][mode] = np.ones(len(setup.mzero[kind][mode]))

            if modes != 'all':
                print(mode, modes)
                if mode.lower() not in modes:
                    continue
            if mode in setup.modes_sc:
                smin = 0
                smax = setup.modes_sc[mode]
            else:
                smin, smax = setup.modes_cc[mode]

            count = 0

            # find correct linear slope
            # xfac = get_m(smin, smax, function)
            print(mode, smin, smax)
            for i in np.arange(smin, smax+2, 2):
                if degs != 'all':
                    if i not in degs:
                        if i == 0:
                            count += 2
                        else:
                            count += 1
                        continue
                if i == 0:
                    for j in [0, 1]:
                        mdamp[kind][mode][count] = get_m(smin, smax,
                                                         function, i)
                        count += 1
                else:
                    for j in np.arange(0, 2*i+1, 1):
                        mdamp[kind][mode][count] = get_m(smin, smax,
                                                         function, i)
                        count += 1

    return mdamp


def calc_Q(mode, fc, ImC00, err=None):
    if err is None:
        _Q = ((mode.freq * 1e3) / (2 * mode.Q)) + \
            1. / np.sqrt(4. * np.pi) * ImC00
    else:
        _Q = ((mode.freq * 1e3) / (2 * mode.Q)) + \
            1. / np.sqrt(4. * np.pi) * (ImC00 + err)
    return 0.5 * fc / _Q
