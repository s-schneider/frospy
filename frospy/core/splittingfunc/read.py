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
from future.utils import native_str
from obspy.core import AttribDict
from obspy.core.event.base import QuantityError

import getpass
import numpy as np
import glob
import os
import subprocess
import pickle
import json

from frospy.util.converter import convert_AB_to_cst
from frospy.util.read import (
    read_modes_in, get_mode_names, get_mode_deg, get_mode_deg_dst, read_pickle
                            )
from frospy.util.base import (chunking_list, split_digit_nondigit,
                            max_cc_degrees, max_sc_degrees, find_unique_name,
                            fQ2cst, fQ2cst_err, get_err)
from frospy import data as frospydata

from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes, Mode, format_name
from frospy.core.database.query import cst_query
from frospy.core.database.write import _write_cst_coeffs


def read(ifile, format=None):
    if format is None:
        # try to guess format from file extension
        _, format = os.path.splitext(ifile)
        format = format[1:]

    if format == 'pickle':
        files = glob.glob(ifile)
        SF = None
        for f in files:
            if SF is None:
                SF = read_pickle(f)
            else:
                SF += read_pickle(f)
        return SF
    return


def read_cst(setup=None, modes=None, cfile=None, modes_dir=None, R=-0.2,
             model='data', include_CRUST=True,
             mdcplbin=None, mdcplccbin=None, verbose=False):
    """
    param modes_dir: path to directory containing:
                         modes.in
                         modes_cc.in
                     as defined for synseis hybrid

    type  modes_dir: string

    Example for new coefficients of hybrid code:def read_cst
    Modes 0S15 - 0T16
    The order of coefficients are as follows:
    indices    mode
    1-231       0 S 15 cst
    232-237     0 T 16 cst
    238-242     0S15-0T16 cst
    243         0 S 15 d00
    244         0 T 16 d00
    """

    if setup is not None or modes_dir is not None:

        if cfile in ('S20RTS', 'S40RTS', 'SP12RTS', 'QRFSI12', 'CRUST',
                     'VSXI', 'VS') or cfile.endswith('.sph'):

            if len(setup.modes_cc) > 0 and len(setup.modes_sc) == 0:
                allmodes = read_modes()
                for _m in setup.modes_cc.keys():
                    _m1, _m2 = _m.split('-')
                    setup.modes_sc[_m1] = 20
                    setup.modes_sc[_m2] = 20
                    setup.modes += allmodes.select(name=_m1)
                    setup.modes += allmodes.select(name=_m2)
        out = read_setup_stats(setup, modes_dir)

        modes_sc, modes_cc, modesin, modes_ccin, modes_scin_dst = out[:]

        sc, cc = get_mode_names(modesin, modes_ccin)
        if cc is not None:
            for m in cc:
                header = {'n': -1, 'type': 'CC', 'l': -1, 'name': m,
                          'sens': None, 'freq': 0, 'Q': 0}
                modes_cc += Mode(header)
    elif modes is not None:
        modes_sc, modes_cc, modesin, modes_ccin = get_modes4cst(modes)
        modes_scin_dst = None
    # cfile in ('S20RTS', 'S40RTS', 'SP12RTS', 'QRFSI12')
    else:
        if not cfile.endswith('sqlite3'):
            print('if cfile not "db", setup or modes_dir has to be given')
            print('trying to load all modes from given cfile')
            modesin = None
            modes_ccin = None
            # return

    if verbose is True:
        print(cfile, modesin, modes_ccin)

    cst, dst, cst_errors, dst_errors = None, None, None, None

    if mdcplbin is not None:
        cst, dst = read_cst_S20RTS(modesin=modesin, modes_ccin=modes_ccin,
                                   setup=setup, modes_dst=modes_scin_dst,
                                   R=R, model=cfile, verbose=verbose,
                                   include_CRUST=include_CRUST,
                                   mdcplbin=mdcplbin,
                                   mdcplccbin=mdcplccbin)
    elif cfile == 'AD':
        file_name = "AD_cst.json"
        path = "%s/AD/%s" % (frospydata.__path__[0], file_name)
        cst, dst, cst_errors, dst_errors = read_cst_AD(modesin, modes_ccin,
                                                       path)
    elif cfile in ['STS_SC', 'STS_GC_SC', 'STS_GC_CC']:
        file_name = "{}.dat".format(cfile)
        path = "%s/STS/%s" % (frospydata.__path__[0], file_name)
        cst, dst, cst_errors, dst_errors = read_cst_AD(modesin, modes_ccin,
                                                       path)

    elif cfile == 'RR':
        cst, dst, cst_errors, dst_errors = read_cst_RR(modesin, modes_ccin,
                                                       verbose=verbose)
    elif cfile in ('S20RTS', 'S40RTS', 'SP12RTS', 'QRFSI12', 'CRUST', 'VSXI', 'VS'):
        cst, dst = read_cst_S20RTS(modesin=modesin, modes_ccin=modes_ccin,
                                   setup=setup, modes_dst=modes_scin_dst,
                                   R=R, model=cfile, verbose=verbose,
                                   include_CRUST=include_CRUST)

    elif cfile == 'TZ':
        cst, dst, cst_errors, dst_errors = read_cst_TZ(modesin, modes_ccin)

    elif cfile == 'REM':
        # temporary solution for REM with measurements in two locations
        cst, cst_errors = read_cst_REM(modesin, modes_ccin)
        if len(cst) == 0:
            try:
                cst, dst, cst_errors, dst_errors = read_cst_db(setup=setup,
                                                               model=cfile)
            except Exception:
                pass
    elif cfile == 'PREM':
        cst, dst = read_cst_PREM(modesin, modes_ccin)

    elif cfile in ['CB', 'TCB']:
        cst = read_cst_RR(modesin, modes_ccin, cfile)

    elif cfile in ['MW', 'WZM']:
        cst, dst, cst_errors, dst_errors = read_cst_MW(modesin, modes_ccin,
                                                       folder_name=cfile)

    elif cfile.upper() == 'SAS':
        cst, dst, cst_errors, dst_errors = read_cst_SAS(modes, setup)

    elif cfile in ['S20RTS+CRUST+BT', 'S20RTS+CRUST+Tr',
                   'S20RTS+CRUST+Ro', 'S20RTS+CRUST+Wh',
                   'BT', 'Tr', 'Ro', 'Wh',
                   'HT', 'QM1', 'DE', 'GLW', 'GD', 'Sumatra',
                   '$\phi$=1.20', '$\phi$=1.10', '$\phi$=1.04',
                   '$\phi$=0.96', '$\phi$=0.90', '$\phi$=0.80']:
        cst, dst, cst_errors, dst_errors = read_cst_db(setup=setup,
                                                       model=cfile)

    elif cfile.endswith('sqlite3'):
        out = read_cst_db(setup=setup, modes=modes, model=model,
                          file_name=cfile)
        cst, dst, cst_errors, dst_errors = out[:]
        modes_sc = Modes()
        modes_cc = Modes()
        if type(modes) is not list:
            modes = [modes]
        for m in modes:
            if '-' not in m:
                modes_sc += read_modes(modenames=m)
            else:
                m_temp = Mode()
                m_temp.name = m
                modes_cc += m_temp
    elif cfile == 'PK':
        cst, dst, cst_errors, dst_errors = read_cst_PK(modesin, modes_ccin)
    else:
        try:
            c, c_tmp, noc = _read_cst_file(cfile, setup)
            # Preparing cst and dst files
            cst, dst = get_cst(modes=modesin, modes_cc=modes_ccin,
                               modes_dst=modes_scin_dst, c=c, noc=noc)

            cst_errors, dst_errors = get_cst_errors(c=c_tmp, modes=modesin,
                                                    modes_cc=modes_ccin,
                                                    modes_dst=modes_scin_dst)
        except ValueError:
            cst, dst, cst_errors, dst_errors = read_cst_AD(modesin, modes_ccin,
                                                           cfile)

    if modesin is None:
        modes_sc = Modes()
        modes_cc = Modes()
        modes_all = read_modes()
        for m in cst.keys():
            modes_sc += modes_all.select(name=m)

    return cst, dst, cst_errors, dst_errors, modes_sc, modes_cc


def _read_cst_file(cfile, setup=None):
    # Read coefficients file
    with open(cfile, 'r') as fh:
        c_tmp = np.genfromtxt(fh)
        if c_tmp.ndim > 1 or setup is not None:
            if c_tmp.ndim > 1:
                # noc is set to False, if Arwens format is read
                noc = False
                c = c_tmp.transpose()[0]
            else:
                # noc is set to False, if Arwens format is read
                noc = False
                c = c_tmp
        else:
            # noc is set to the number of cst, if Haydars format is read
            noc = int(c_tmp[0])
            c = c_tmp[1:]
    return c, c_tmp, noc


def read_setup_stats(setup, modes_dir):
    # Read mode files
    modes_scin_dst = None
    if setup is not None:
        modesin, modes_ccin = setup.get_modes_in()
        if hasattr(setup, 'get_modes_in_dst'):
            if setup.get_modes_in_dst() is not None:
                modes_scin_dst = setup.get_modes_in_dst()
    elif modes_dir is not None:
        modesin, modes_ccin = read_modes_in(modes_dir)
    modes_sc = Modes()
    modes_cc = Modes()
    modes_all = read_modes()
    sc, cc = get_mode_names(modesin, modes_ccin)
    for m in sc:
        modes_sc += modes_all.select(name=m)

    return (modes_sc, modes_cc, modesin, modes_ccin,
            modes_scin_dst)


def get_cst_errors(c, modes, modes_cc, modes_dst=None, modes_cc_dst=None):
    # to check size of mcst
    if not isinstance(c, list):
        try:
            if c.shape[1] == 1:
                c = np.array([c])
            else:
                c = np.array([x for x in zip(c)])  # x = zip(*c)
        except (IndexError, AttributeError):
            return None, None

        if c.shape[1] == 2 or c.shape[1] == 3:  # standard input
            # c.ndim == 2 # standard input
            # c.ndim == 3 # Arwens format final format
            c = c.transpose()[-1]
        else:
            # c.ndim == 1 # simple cst input
            # c.ndim == 4 # current inversion output
            cst_errors = None
            dst_errors = None
            return cst_errors, dst_errors
    ind = 0
    count = 0
    cst_errors = AttribDict()
    dst_errors = AttribDict()

    for mode in modes[1:]:
        count += 1
        m = mode.split()
        name = ''.join(m[0:3])
        max_cdeg = int(m[3])
        # max_ddeg = int(m[4])
        cdegs = np.arange(2, max_cdeg+2, 2)

        # only read if mode has smax > 0 and radial modes
        if ((max_cdeg > 0 and int(m[2]) > 0) or (m[2] == '0')):
            cst_errors[name] = AttribDict()
            dst_errors[name] = AttribDict()
            # c00 coefficients
            err = QuantityError()
            err.uncertainty = np.array([c[ind]])
            err.upper_uncertainty = np.array([c[ind]])
            err.lower_uncertainty = np.array([c[ind]])
            err.confidence_level = 0
            cst_errors[name]['0'] = err
            ind += 1

            # d00 coefficients
            err = QuantityError()
            err.uncertainty = np.array([c[ind]])
            err.upper_uncertainty = np.array([c[ind]])
            err.lower_uncertainty = np.array([c[ind]])
            err.confidence_level = 0
            dst_errors[name]['0'] = err
            ind += 1
            for d in cdegs:
                ind_end = ind + 2*d + 1
                err = QuantityError()
                err.uncertainty = np.array(c[ind:ind_end])
                err.upper_uncertainty = np.array(c[ind:ind_end])
                err.lower_uncertainty = np.array(c[ind:ind_end])
                err.confidence_level = 0
                cst_errors[name][str(d)] = err
                ind = ind_end

        # Cross coupling coefficients
        # CC after each mode in Self
        # 1-1  ->  1-2  ->  1-3
        #          2-2  ->  2-3
        #                   3-3
        if modes_cc is not None:
            # reading for 2 and 3 modes
            if count == 1 and len(modes_cc) >= 2:
                # 1-2  ->  1-3
                cc = modes_cc[1:3]
            elif count == 2 and len(modes_cc) == 4:
                # ->  2-3
                cc = modes_cc[-1:]
            else:
                cc = []

            for mode_cc in cc:
                m = mode_cc.split()
                namecc = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                ccdegs = max_cc_degrees(m[:-2])
                max_cdeg = int(m[-1])
                min_cdeg = min(ccdegs)

                if min_cdeg != int(m[-2]) and max_cdeg > 0:
                    i = list(ccdegs).index(int(m[-2]))
                    ccdegs = ccdegs[i:]

                if max_cdeg > 0:  # only read if CC has smax > 0
                    cst_errors[namecc] = AttribDict()
                    dst_errors[namecc] = AttribDict()
                    for d in ccdegs:
                        if d > max_cdeg:
                            break
                        if d == 0:
                            # c00 coefficients
                            # instead of coefss use the Error class here
                            err = QuantityError()
                            err.uncertainty = np.array(c[ind])
                            err.upper_uncertainty = np.array(c[ind])
                            err.lower_uncertainty = np.array(c[ind])
                            err.confidence_level = 0
                            cst_errors[namecc]['0'] = err
                            ind += 1

                            # d00 coefficients
                            err = QuantityError()
                            err.uncertainty = np.array(c[ind])
                            err.upper_uncertainty = np.array(c[ind])
                            err.lower_uncertainty = np.array(c[ind])
                            err.confidence_level = 0
                            dst_errors[namecc]['0'] = err
                            ind += 1
                        else:
                            ind_end = ind + 2*d + 1
                            err = QuantityError()
                            err.uncertainty = np.array(c[ind:ind_end])
                            err.upper_uncertainty = np.array(c[ind:ind_end])
                            err.lower_uncertainty = np.array(c[ind:ind_end])
                            err.confidence_level = 0
                            cst_errors[namecc][str(d)] = err
                            ind = ind_end
    # dst coeffs
    count = 0
    if modes_dst is not None:
        # loop through dst
        for mode in modes_dst[1:]:
            count += 1
            m = mode.split()
            name = ''.join(m[0:3])
            max_cdeg = int(m[3])
            cdegs = np.arange(2, max_cdeg+2, 2)

            # only read if mode has smax > 2
            if max_cdeg >= 2:
                for d in cdegs:
                    ind_end = ind + 2*d + 1
                    err = QuantityError()
                    err.uncertainty = np.array(c[ind:ind_end])
                    err.upper_uncertainty = np.array(c[ind:ind_end])
                    err.lower_uncertainty = np.array(c[ind:ind_end])
                    err.confidence_level = 0
                    dst_errors[name][str(d)] = err
                    ind = ind_end

            # Cross coupling coefficients
            # CC after each mode in Self
            # 1-1  ->  1-2  ->  1-3
            #          2-2  ->  2-3
            #                   3-3
            if modes_cc_dst is not None:
                # reading for 2 and 3 modes
                if count == 1 and len(modes_cc_dst) >= 2:
                    # 1-2  ->  1-3
                    cc = modes_cc_dst[1:3]
                elif count == 2 and len(modes_cc_dst) == 4:
                    # ->  2-3
                    cc = modes_cc_dst[-1:]
                else:
                    cc = []

                for mode_cc in cc:
                    m = mode_cc.split()
                    namecc = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                    ccdegs = max_cc_degrees(m[:-2])
                    max_cdeg = int(m[-1])
                    min_cdeg = min(ccdegs)

                    if min_cdeg != int(m[-2]) and max_cdeg > 0:
                        i = list(ccdegs).index(int(m[-2]))
                        ccdegs = ccdegs[i:]

                    if max_cdeg > 0:  # only read if CC has smax > 0
                        for d in ccdegs:
                            if d > max_cdeg:
                                break

                            ind_end = ind + 2*d + 1
                            err = QuantityError()
                            err.uncertainty = np.array(c[ind:ind_end])
                            err.upper_uncertainty = np.array(c[ind:ind_end])
                            err.lower_uncertainty = np.array(c[ind:ind_end])
                            err.confidence_level = 0
                            dst_errors[namecc][str(d)] = err
                            ind = ind_end

    return cst_errors, dst_errors


def get_cst(modes, modes_cc, c, noc, modes_dst=None, modes_cc_dst=None):
    # Arwens format
    if noc is False:
        ind = 0
        count = 0
        cst = AttribDict()
        dst = AttribDict()

        for mode in modes[1:]:
            count += 1
            m = mode.split()
            name = ''.join(m[0:3])
            max_cdeg = int(m[3])
            # max_ddeg = int(m[4])
            cdegs = np.arange(2, max_cdeg+2, 2)

            # only read if mode has smax > 0 and radial modes
            if ((max_cdeg > 0 and int(m[2]) > 0) or (m[2] == '0')):
                cst[name] = AttribDict()
                dst[name] = AttribDict()
                # c00 coefficients
                cst[name]['0'] = np.array([c[ind]])
                ind += 1

                # d00 coefficients
                dst[name]['0'] = np.array([c[ind]])
                ind += 1
                for d in cdegs:
                    ind_end = ind + 2*d + 1
                    # Check here if file has actually values for that deg
                    # Only for data/PK files
                    if ind_end > len(c):
                        break
                    cst[name][str(d)] = np.array(c[ind:ind_end])
                    ind = ind_end

            # Cross coupling coefficients
            # CC after each mode in Self
            # 1-1  ->  1-2  ->  1-3
            #          2-2  ->  2-3
            #                   3-3
            if modes_cc is not None:
                # reading for 2 and 3 modes
                if count == 1 and len(modes_cc) >= 2:
                    # 1-2  ->  1-3
                    cc = modes_cc[1:3]
                elif count == 2 and len(modes_cc) == 4:
                    # ->  2-3
                    cc = modes_cc[-1:]
                else:
                    cc = []

                for mode_cc in cc:
                    m = mode_cc.split()
                    namecc = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                    ccdegs = max_cc_degrees(m[:-2])
                    max_cdeg = int(m[-1])
                    min_cdeg = min(ccdegs)

                    if min_cdeg != int(m[-2]) and max_cdeg > 0:
                        i = list(ccdegs).index(int(m[-2]))
                        ccdegs = ccdegs[i:]

                    if max_cdeg > 0:  # only read if CC has smax > 0
                        cst[namecc] = AttribDict()
                        dst[namecc] = AttribDict()
                        for d in ccdegs:
                            if d > max_cdeg:
                                break
                            if d == 0:
                                # c00 coefficients
                                cst[namecc]['0'] = np.array([c[ind]])
                                ind += 1

                                # d00 coefficients
                                dst[namecc]['0'] = np.array([c[ind]])
                                ind += 1
                            else:
                                ind_end = ind + 2*d + 1
                                cst[namecc][str(d)] = np.array(c[ind:ind_end])
                                ind = ind_end
        # dst coeffs
        count = 0
        if modes_dst is not None:
            # loop through dst
            for mode in modes_dst[1:]:
                count += 1
                m = mode.split()
                name = ''.join(m[0:3])
                max_cdeg = int(m[3])
                cdegs = np.arange(2, max_cdeg+2, 2)

                # only read if mode has smax > 2
                if max_cdeg >= 2:
                    for d in cdegs:
                        ind_end = ind + 2*d + 1
                        dst[name][str(d)] = np.array(c[ind:ind_end])
                        ind = ind_end

                # Cross coupling coefficients
                # CC after each mode in Self
                # 1-1  ->  1-2  ->  1-3
                #          2-2  ->  2-3
                #                   3-3
                if modes_cc_dst is not None:
                    # reading for 2 and 3 modes
                    if count == 1 and len(modes_cc_dst) >= 2:
                        # 1-2  ->  1-3
                        cc = modes_cc_dst[1:3]
                    elif count == 2 and len(modes_cc_dst) == 4:
                        # ->  2-3
                        cc = modes_cc_dst[-1:]
                    else:
                        cc = []

                    for mode_cc in cc:
                        m = mode_cc.split()
                        namecc = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                        ccdegs = max_cc_degrees(m[:-2])
                        max_cdeg = int(m[-1])
                        min_cdeg = min(ccdegs)

                        if min_cdeg != int(m[-2]) and max_cdeg > 0:
                            i = list(ccdegs).index(int(m[-2]))
                            ccdegs = ccdegs[i:]

                        if max_cdeg > 0:  # only read if CC has smax > 0
                            for d in ccdegs:
                                if d > max_cdeg:
                                    break
                                ind_end = ind + 2*d + 1
                                dst[namecc][str(d)] = np.array(c[ind:ind_end])
                                ind = ind_end
    else:
        # Haydars format
        # Self coupling coefficients, dst are commented out
        ind = 0
        cst = AttribDict()
        for mode in modes[1:]:
            m = mode.split()
            name = ''.join(m[0:3])
            cst[name] = AttribDict()

            max_cdeg = int(m[3])
            cdegs = max_sc_degrees(int(m[2]))

            for d in cdegs:
                if d > max_cdeg:
                    break
                ind_end = ind + 2*d + 1
                cst[name][str(d)] = np.array(c[ind:ind_end])
                ind = ind_end

        if modes_cc is not None:
            # Cross coupling coefficients
            for mode in modes_cc[1:]:
                m = mode.split()
                name = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                cst[name] = AttribDict()

                ccdegs = max_cc_degrees(m[:-2])
                max_cdeg = int(m[-2])

                for d in ccdegs:
                    if d > max_cdeg:
                        break
                    ind_end = ind + 2*d + 1
                    cst[name][str(d)] = np.array(c[ind:ind_end])
                    ind = ind_end

        dst = AttribDict()
        for mode in modes[1:]:
            m = mode.split()
            name = ''.join(m[0:3])
            dst[name] = AttribDict()

            max_ddeg = int(m[4])
            ddegs = max_sc_degrees(int(m[2]))

            for d in ddegs:
                if d > max_ddeg:
                    break
                ind_end = ind + 2*d + 1
                dst[name][str(d)] = np.array(c[ind:ind_end])
                ind = ind_end

        if modes_cc is not None:
            # Cross coupling coefficients
            for mode in modes_cc[1:]:
                m = mode.split()
                name = ''.join(m[0:3]) + '-' + ''.join(m[3:6])
                dst[name] = AttribDict()

                ccdegs = max_cc_degrees(m[:-2])
                max_ddeg = int(m[-1])

                for d in ccdegs:
                    if d > max_ddeg:
                        break
                    ind_end = ind + 2*d + 1
                    dst[name][str(d)] = np.array(c[ind:ind_end])
                    ind = ind_end
                if len(dst[name]) == 0:
                    dst.pop(name, None)

    return cst, dst


def get_cst_dat(meas, cst, dst, cst_errors, dst_errors, mnames):
    mode = meas[0].split()
    c_type = 'sc'
    if len(mode) == 2:
        _m = format_name("{}S{}".format(mode[0], mode[1]))
    elif len(mode) == 3:
        _m = format_name("{}{}{}".format(mode[0], mode[1], mode[2]))
    elif len(mode) == 6:
        _m = format_name("{}{}{}".format(mode[0], mode[1], mode[2]))
        _m += '-'
        _m += format_name("{}{}{}".format(mode[3], mode[4], mode[5]))
        c_type = 'cc'
    if mnames is not None:
        if _m not in mnames[0]:
            if _m not in mnames[1]:
                return cst, dst, cst_errors, dst_errors
    if _m not in cst:
        cst[_m] = {}
        cst_errors[_m] = {}
        dst[_m] = {}
        dst_errors[_m] = {}

    if c_type == 'sc':
        coeffs = chunking_list(meas[1:], 2)
    else:
        coeffs = chunking_list(meas[2:], 2)
    for c in coeffs:
        _cst = c[0].split()
        _err = c[1].split()

        if len(_cst) == 2:
            deg = '0'
            cst[_m][deg] = np.array([float(_cst[0])])
            cst_errors[_m][deg] = get_err(float(_err[0]))
            dst[_m][deg] = np.array([float(_cst[1])])
            dst_errors[_m][deg] = get_err(float(_err[1]))
        else:
            deg = str(int((len(_cst) - 1) / 2))
            x = [float(y) for y in _cst]
            cst[_m][deg] = np.array(x)
            x = [float(y) for y in _cst]
            cst_errors[_m][deg] = get_err(x)

    return cst, dst, cst_errors, dst_errors


def read_cst_SAS(modes, setup):
    db_file = "%s/SAS/cst.sqlite3" % frospydata.__path__[0]
    cst, dst, cst_errors, dst_errors = read_cst_db(setup=setup, modes=modes,
                                                   file_name=db_file,
                                                   model='data')
    return cst, dst, cst_errors, dst_errors


def read_cst_db(model, setup=None, modes=None, file_name=None, R=-0.2):
    """
    read cst and dst coeffs from database.
    Default model is 'data' for measured cst
    """
    if file_name is None:
        file_name = "%s/cst.sqlite3" % frospydata.__path__[0]

    if setup is None and modes is None:
        # msg = 'No input specified'
        # print(msg)
        return

    coeffs = AttribDict({'cst': {}, 'dst': {}})
    coeffs_err = AttribDict({'cst': {}, 'dst': {}})

    # modes in sc
    if setup is not None:
        m = setup.modes_sc.keys()
        modes = [format_name(x).upper() for x in m]
        # modes in cc
        m = setup.modes_cc.keys()
        if len(m) > 0:
            modes_cc = [format_name(x).upper() for x in m]
            modes.extend(modes_cc)

    dbq = cst_query(db_path=file_name, model=model, modes=modes)

    for line in dbq:
        err = QuantityError()
        mode = str(line[0])
        kind = str(line[1])
        deg = line[2]
        x = [float(y.split()[0]) for y in line[3].split(',')]
        xe = [float(y.split()[1]) for y in line[3].split(',')]

        try:
            xe_u = [float(y.split()[2]) for y in line[3].split(',')]
            xe_l = [float(y.split()[3]) for y in line[3].split(',')]
            xe_c = [float(y.split()[4]) for y in line[3].split(',')]
            err.upper_uncertainty = np.array(xe_u)
            err.lower_uncertainty = np.array(xe_l)
            err.confidence_level = xe_c
        except IndexError:
            pass

        err.uncertainty = np.array(xe)

        if mode not in coeffs[kind]:
            coeffs[kind][mode] = AttribDict()
            coeffs_err[kind][mode] = AttribDict()
        coeffs[kind][mode][deg] = np.array(x)
        coeffs_err[kind][mode][deg] = err

    # Check here for missing modes and degrees, if model is S20RTS
    #     if mode not in modes_checked:
    #         modes_checked[mode] = int(deg)
    #
    # nodb_modes = list(set(modes) - set(modes_checked))

    # Write missing modes into db
    # if model in ['S20RTS']:
        # Calculate missing modes and degrees
        # cst, dst = read_cst_S20RTS(modesin=modesin, modes_ccin=modes_ccin,
        #                            setup=setup, modes_dst=modes_scin_dst,
        #                            R=R)
        # for mode in nodb_modes:
        #     continue

    return coeffs['cst'], coeffs['dst'], coeffs_err['cst'], coeffs_err['dst']


def read_cst_MW(modesin, modes_ccin, folder_name="MW"):

    path = "%s/%s/fc_and_Q.txt" % (frospydata.__path__[0], folder_name)

    cst = AttribDict()
    cst_errors = AttribDict()
    dst = AttribDict()
    dst_errors = AttribDict()
    allmodes = read_modes()

    mnames = get_mode_names(modesin, modes_ccin)[0]
    # modes_prem = read_mode_class()
    with open(path) as fh:
        content = fh.readlines()

    for line in content:
        if line.startswith('#'):
            continue
        name = ''.join( line.split()[:3])
        if name not in mnames:
            continue
        mode = allmodes.select(name=name)[0]

        # w_prem = modes_prem.select(name=name)[0].freq
        f0 = float(line.split()[3])
        f0_err = float(line.split()[4])

        if name not in cst:
            cst[name] = AttribDict()
            cst_errors[name] = AttribDict()
        if name not in dst:
            dst[name] = AttribDict()
            dst_errors[name] = AttribDict()

        c00 = fQ2cst(f0, mode.Q, mode)
        cst[name]['0'] = np.array([c00[0]],
                                  dtype='float')
        dst[name]['0'] = np.array([c00[1]],
                                  dtype='float')

        c00 = fQ2cst_err(f0, f0_err, mode.Q, 0, mode,
                         c00[0], c00[1])
        cst_errors[name]['0'] = get_err([c00[0]])
        dst_errors[name]['0'] = get_err([c00[1]])

    return cst, dst, cst_errors, dst_errors


def read_cst_AD(modesin, modes_ccin, file_name):
    if file_name.endswith('json'):
        return read_cst_AD_json(modesin, modes_ccin, file_name)
    else:
        cst = AttribDict()
        cst_errors = AttribDict()
        dst = AttribDict()
        dst_errors = AttribDict()
        if modesin:
            mnames = get_mode_names(modesin, modes_ccin)
        else:
            mnames = None

        with open(file_name, 'r') as fh:
            content = fh.readlines()

        # Skip first i lines, they are comments
        for i, line in enumerate(content):
            if line[0].isnumeric():
                istart = i
                break

        meas = []
        for line in content[istart:]:
            if line[0].isnumeric() and not line[2].isnumeric():
                if len(meas) != 0:
                    cst, dst, cst_errors, dst_errors = get_cst_dat(meas, cst,
                                                                   dst,
                                                                   cst_errors,
                                                                   dst_errors,
                                                                   mnames)
                meas = []
            meas.append(line)

        return cst, dst, cst_errors, dst_errors


def read_cst_AD_json(modes, modes_cc, file_name="AD_cst.json"):
    # path = "%s/AD/%s" % (frospydata.__path__[0], file_name)


    with open(file_name, 'r') as fh:
        data = json.load(fh)

    mnames = get_mode_names(modes, modes_cc)
    cst = None
    dst = None
    cst_errors = None
    dst_errors = None
    for mode, values in data['SF'].items():
        name = format_name(mode.upper())
        if name not in mnames[0]:
            continue

        lcut = values['lcut']
        _n = split_digit_nondigit(name)
        _m = ['1', "%i %s %i %i %i" % (int(_n[0]), _n[1], int(_n[2]), lcut, 0)]
        cst, dst = get_cst(modes=_m, modes_cc=None, c=values['cst'], noc=False)
        cst_errors, dst_errors = get_cst_errors(c=values['cst_errors'],
                                                modes=_m,
                                                modes_cc=None)
    return cst, dst, cst_errors, dst_errors


def read_cst_PREM(modesin, modes_ccin):

    sc_modes, cc_modes = get_mode_names(modesin, modes_ccin)
    sc_cdeg, sc_ddeg, cc_cdeg, cc_ddeg = get_mode_deg(modesin, modes_ccin)

    cst = AttribDict()
    for mode, s_max in zip(sc_modes, sc_cdeg):
        for deg in range(0, int(s_max) + 1, 2):
            if mode not in cst:
                cst[mode] = AttribDict()
            cst[mode][str(deg)] = np.zeros(2 * deg + 1)
    dst = None
    return cst, dst


def read_cst_TZ(modes, modes_cc, file_name="AB_s2_tromp_zanzerkia.dat",
                sdegs=['0', '2']):
    """
    Possible filenames:
            AB_s2_tromp_zanzerkia.dat
    """

    data_list = ['AB_s2_tromp_zanzerkia.dat']
    if file_name in data_list:
        path = "%s/TZ/%s" % (frospydata.__path__[0], file_name)
    else:
        path = file_name
    cst = AttribDict()
    cst_errors = AttribDict()
    dst = AttribDict()
    dst_errors = AttribDict()
    allmodes = read_modes()

    mnames = get_mode_names(modes, modes_cc)[0]
    # modes_prem = read_mode_class()
    with open(path) as fh:
        content = fh.readlines()

    # first line in file is the unit of the coefficients
    exp = float(content[0])
    content = chunking_list(content[1:], 4)
    for line in content:
        name = line[0].rstrip().upper()
        if name not in mnames:
            continue
        mode = allmodes.select(name=name)[0]

        # w_prem = modes_prem.select(name=name)[0].freq
        f0 = float(line[1].split()[0])
        Q = float(line[1].split()[1])
        f0_err = float(line[2].split()[0])
        Q_err = float(line[2].split()[1])

        L = np.array(line[1].split()[1:]).astype(float)
        L_err = np.array(line[2].split()[1:]).astype(float)
        # 1E-6 is the unit presented in the paper
        A = np.array([L[1], L[2], L[4]]) * exp
        B = np.array([L[3], L[5]]) * exp

        A_err = np.array([L_err[1], L_err[2], L_err[4]]) * exp
        B_err = np.array([L_err[3], L_err[5]]) * exp

        if name not in cst:
            cst[name] = AttribDict()
            cst_errors[name] = AttribDict()
        if name not in dst:
            dst[name] = AttribDict()
            dst_errors[name] = AttribDict()

        for deg in sdegs:
            if deg == '0':
                c00 = fQ2cst(f0, Q, mode)
                cst[name][deg] = np.array([c00[0]],
                                          dtype='float')
                dst[name][deg] = np.array([c00[1]],
                                          dtype='float')

                c00 = fQ2cst_err(f0, f0_err, Q, Q_err, mode,
                                 c00[0], c00[1])
                cst_errors[name][deg] = get_err([c00[0]])
                dst_errors[name][deg] = get_err([c00[1]])
                continue
            void, mcst = convert_AB_to_cst(A, B, int(deg))
            cst[name][deg] = mcst * f0

            void, mcst_err = convert_AB_to_cst(A + A_err, B + B_err, int(deg))
            cst_errors[name][deg] = get_err(abs(mcst - mcst_err) * f0)
            dst_errors[name][deg] = get_err(abs(mcst - mcst_err) * f0,
                                            zeros=True)
    return cst, dst, cst_errors, dst_errors


def read_cst_PK(modesin, modes_ccin, verbose=False):
    """
    Comment for get_cst in this function:
    Default values modes_ccin is None, because no CC values are in
    our database
    """
    path = frospydata.__path__[0]
    mnames = get_mode_names(modesin, modes_ccin)
    cst = AttribDict()
    cst_errors = AttribDict()
    dst = AttribDict()
    dst_errors = AttribDict()
    for i, mode in enumerate(mnames[0]):
        cfile = "%s/PK/%s.cst" % (path, format_name(mode, 2).lower())
        if not os.path.exists(cfile):
            continue

        c, c_tmp, noc = _read_cst_file(cfile, None)
        if verbose:
            print(mode)

        _cst, _dst = get_cst(modes=[1, modesin[i + 1]], modes_cc=None, c=c,
                             noc=False, modes_dst=None)
        cst[mode], dst[mode] = _cst[mode], _dst[mode]

        _cst_err = c_tmp.transpose()[2]
        _cst_err, _dst_err = get_cst(modes=[1, modesin[i + 1]], modes_cc=None,
                                     c=_cst_err, noc=False, modes_dst=None)

        for deg, _c in _cst_err[mode].items():
            if mode not in cst_errors:
                cst_errors[mode] = AttribDict()
            if mode not in dst_errors:
                dst_errors[mode] = AttribDict()
            cst_errors[mode][deg] = get_err(_c)
            dst_errors[mode][deg] = get_err(_c, zeros=True)

    return cst, dst, cst_errors, dst_errors


def read_cst_RR(modesin, modes_ccin, mname='RR', verbose=False):
    """
    Docstring will follow
    """
    # This is the reading function
    def read_cst_errors_RR(content, cst, cst_errors, dst, dst_errors,
                           fQ, errors, mnames, deg, c_type, allmodes,
                           verbose):
        if c_type == 'sc':
            for line in content:
                n = ''.join(line.split()[:3])
                name = format_name(n.upper())

                if verbose is True:
                    print('Searching SC: ', name, mnames[0])
                if format_name(name) not in mnames[0]:
                    continue
                if verbose is True:
                    print('Found: ', name)
                mode = allmodes.select(name=name)[0]
                if name not in cst:
                    cst[name] = AttribDict()
                    cst_errors[name] = AttribDict()
                if name not in dst:
                    dst[name] = AttribDict()
                    dst_errors[name] = AttribDict()
                if deg == '0':
                    f0 = float(line.split()[3])
                    Q = float(line.split()[4])
                    fQ[name] = [f0, Q]
                    c00 = fQ2cst(f0, Q, mode)
                    cst[name][deg] = np.array([c00[0]],
                                              dtype='float')
                    cst_errors[name][deg] = get_err([0], zeros=True)
                    dst[name][deg] = np.array([c00[1]],
                                              dtype='float')
                    dst_errors[name][deg] = get_err([0], zeros=True)
                else:
                    if errors is False:
                        cst[name][deg] = np.array(line.split()[3:],
                                                  dtype='float')
                    else:
                        val = np.array(line.split()[3:], dtype='float')
                        cst_errors[name][deg] = get_err(val)

        elif c_type == 'cc' and mnames[1] is not None:
            for line in content:
                name1 = ''.join(line.split()[:3]).upper()
                name2 = ''.join(line.split()[3:6]).upper()
                name = '-'.join([name1, name2])
                name_r = '-'.join([name2, name1])

                mt1 = split_digit_nondigit(name1)[1]
                mt2 = split_digit_nondigit(name2)[1]

                if mt1.upper() != mt2.upper():
                    ccST = True
                else:
                    ccST = False

                if verbose is True:
                    print('Searching CC: ', name, mnames[1])

                if name not in mnames[1] and name_r not in mnames[1]:
                    continue

                if verbose is True:
                    print('Found: ', name, 'ST crosscoupling: ', ccST)

                if ccST is True:
                    ccfac = -1.
                else:
                    ccfac = 1.
                if name_r in mnames[1]:
                    name = name_r

                if name not in cst:
                    cst[name] = AttribDict()
                    cst_errors[name] = AttribDict()
                if errors is False:
                    cst[name][deg] = ccfac * np.array(line.split()[6:],
                                                      dtype='float')
                else:
                    val = ccfac * np.array(line.split()[6:], dtype='float')
                    cst_errors[name][deg] = get_err(val)
                # cst[name][deg] = np.array(line.split()[6:], dtype='float')
        return cst, cst_errors, dst, dst_errors

    path = frospydata.__path__[0]
    # sc_files = glob.glob('%s/RR*sc.cst' % path)
    # cc_files = glob.glob('%s/RR*cc.cst' % path)
    files = glob.glob('%s/%s/%s*.cst' % (path, mname, mname))
    cst = AttribDict()
    dst = AttribDict()
    fQ = AttribDict()
    cst_errors = AttribDict()
    dst_errors = AttribDict()
    allmodes = read_modes()
    mnames = get_mode_names(modesin, modes_ccin)

    if verbose is True:
        print(mnames)
    # from IPython import embed; embed()
    for f in files:
        if verbose is True:
            print(f)
        with open(f, 'r') as fh:
            content = fh.readlines()

        if mname == 'RR':
            deg = f.split('/')[-1].split('RR')[1].split('-')[0]
            c_type = f.split('/')[-1].split('-')[1].split('.')[0]
        elif mname in ['CB', 'TCB']:
            deg = '2'
            c_type = 'cc'
        if verbose is True:
            print(deg, c_type)
        cst, cst_errors, dst, dst_errors = read_cst_errors_RR(content, cst,
                                                              cst_errors,
                                                              dst,
                                                              dst_errors,
                                                              fQ,
                                                              False,
                                                              mnames,
                                                              deg,
                                                              c_type,
                                                              allmodes,
                                                              verbose)
    # Reading fcenter errorbars:
    if mname == 'RR':
        with open('%s/RR/fQ.dat' % path, 'r') as fh:
            content = fh.readlines()
        for line in content:
            if line.startswith('#'):
                continue
            name1 = ''.join(line.split()[0]).upper()
            mtype = ''.join(line.split()[1]).upper()
            name2 = ''.join(line.split()[2]).upper()
            name = format_name(''.join([name1, mtype, name2]))
            mode = allmodes.select(name=name)[0]
            if name not in mnames[0]:
                continue
            f0_err = float(line.split()[4])
            Q_err = float(line.split()[6])
            f = fQ[name][0]
            Q = fQ[name][1]
            c00 = fQ2cst_err(f, f0_err, Q, Q_err, mode,
                             cst[name]['0'], dst[name]['0'])

            cst_errors[name] = AttribDict()
            dst_errors[name] = AttribDict()
            # c00 coefficients
            cst_errors[name]['0'] = get_err(c00[0])
            dst_errors[name]['0'] = get_err(c00[1])
            for degs, val in cst[name].items():
                if degs == '0':
                    continue
                cst_errors[name][degs] = get_err(val, zeros=True)
                dst_errors[name][degs] = get_err(val, zeros=True)

        # Read the remaining errorbars from file
        files = glob.glob('%s/%s/*err' % (path, mname))
        for f in files:
            if verbose is True:
                print(f)
            with open(f, 'r') as fh:
                content = fh.readlines()

            if mname == 'RR':
                base = os.path.basename(f).split('d')[1]
                deg = base.split('.')[0]
                c_type = base.split('.')[1].split('err')[0]
            elif mname in ['CB', 'TCB']:
                deg = '2'
                c_type = 'cc'
            if verbose is True:
                print(deg, c_type)

            out = read_cst_errors_RR(content, cst, cst_errors,
                                     dst, dst_errors, fQ, True, mnames,
                                     deg, c_type, allmodes, verbose)
            cst, cst_errors, dst, dst_errors = out[:]
    return cst, dst, cst_errors, dst_errors


def read_cst_REM(modesin, modes_ccin):

    def read_REM_lines(content, count, mnames):
        n, l, no_s = np.array(content[count].split(), dtype=int)
        name = "%s%s%s" % (n, mtype, l)
        if name not in mnames:
            count += no_s * 2 + 1
            return cst, cst_errors, count
        if name not in cst:
            cst[name] = AttribDict()
            cst_errors[name] = AttribDict()
        count += 1
        for _i in range(no_s):
            modedata = content[count]
            errordata = content[count + 1]
            deg = modedata.split()[0]
            cst[name][deg] = np.array(modedata.split()[1:], dtype=float)

            err = QuantityError()
            err.uncertainty = np.array(errordata.split(), dtype=float)
            err.upper_uncertainty = np.array(errordata.split(), dtype=float)
            err.lower_uncertainty = np.array(errordata.split(), dtype=float)
            err.confidence_level = 0
            cst_errors[name][deg] = err
            count += 2
        return cst, cst_errors, count

    path = frospydata.__path__[0]
    files = glob.glob('%s/REM/*.txt' % path)
    cst = AttribDict()
    cst_errors = AttribDict()

    mnames = get_mode_names(modesin, modes_ccin)[0]

    for f in files:
        if 'bestcofft.txt' in f:
            mtype = 'T'
        else:
            mtype = 'S'
        with open(f, 'r') as fh:
            content = fh.readlines()

        count = 0
        while True:
            # Read as long as the are lines in file
            try:
                cst, cst_errors, count = read_REM_lines(content, count, mnames)
            except IndexError:
                pass
    # files = glob.glob('%s/REM/fQ?.dat' % path)
    # for f in files:
    # cst, dst, cst_errors, dst_errors = read_REM_fQ()
    return cst, cst_errors


def _read_cst_S20RTS_db(setup, file_name="S20RTS_CRUST.sqlite3",
                        verbose=False):
    bins = ['/quanta1/home', '/net/home']
    for path in bins:
        if os.path.exists(path):
            bin_path = path
            break

    if bin_path is not None:
        if bin_path.startswith('/quanta'):
            path = "/quanta1/home/simons/dev/python/frospy/frospy"
        else:
            if getpass.getuser() == 'simons':
                path = "/net/home/simons/dev/python/frospy/frospy"
            else:
                path = "/net/home/talavera/codes/nmPy/nmpy"
        path = os.path.join(path, 'data')
    else:
        path = frospydata.__path__[0]
    if verbose is True:
        print('read', path, file_name)
    if file_name == "S20RTS_CRUST.sqlite3":
        path = "%s/S20RTS/%s" % (path, file_name)
        out = read_cst_db(setup=setup, model='S20RTS', file_name=path)
    elif file_name == "S40RTS_CRUST.sqlite3":
        path = "%s/S40RTS/%s" % (path, file_name)
        out = read_cst_db(setup=setup, model='S40RTS', file_name=path)
    cst, dst, cst_errors, dst_errors = out[:]

    return cst, dst


def _write_cst_S20RTS_db(cst, dst, file_name="S20RTS_CRUST.sqlite3",
                         verbose=False):
    bins = ['/quanta1/home', '/net/home']
    for path in bins:
        if os.path.exists(path):
            bin_path = path
            break
    if bin_path is None:
        raise IOError('Bin_path not valid')

    if bin_path.startswith('/quanta'):
        path = "/quanta1/home/simons/dev/python/frospy/frospy"
    else:
        if getpass.getuser() == 'simons':
            path = "/net/home/simons/dev/python/frospy/frospy"
        else:
            path = "/net/home/talavera/codes/nmPy/nmpy"

    if file_name == "S20RTS_CRUST.sqlite3":
        path = "%s/data/S20RTS/%s" % (path, file_name)
        model = 'S20RTS'
    elif file_name == "S40RTS_CRUST.sqlite3":
        path = "%s/data/S40RTS/%s" % (path, file_name)
        model = 'S40RTS'
    else:
        model = file_name.split('.')[0]
        path = "{}/data/{}/{}".format(path, model, file_name)
    if verbose is True:
        print('write', path)
    _write_cst_coeffs(cst, dst, path, model=model, author=None, lcut='all')

    return


def read_cst_S20RTS(modesin, modes_ccin, setup=None, bin_path=None,
                    keep_mcst=False, modes_dst=None, modes_cc_dst=None,
                    R=-0.2, model='S20RTS', include_CRUST=True,
                    verbose=False,
                    mdcplbin=None,
                    mdcplccbin=None):
    """
    Calculates S20RTS coefficients using the program defined in
    S20RTS_path for given self-coupling modes in 'modes' and cross-coupling
    modes in 'modes_cc'

    8.6.18: Currently only self-coupling
    """
    # Try to read coefficients from database, if mode is not yet in there or
    # higher degree than in db is requested it will be calculated
    # calculate R dst predictions if not R=-2. Only R=-0.2 saved in database
    # print(model, 'CRUST', include_CRUST)

    # from IPython import embed; embed()
    if R == -0.2:
        try:
            if model == 'S20RTS':
                file_name = "S20RTS_CRUST.sqlite3"
            elif model == 'S40RTS':
                file_name = "S40RTS_CRUST.sqlite3"
            else:
                file_name = "{}.sqlite3".format(model)
            cst, dst = _read_cst_S20RTS_db(setup, file_name=file_name)
            if len(cst) == 0:
                raise IOError

            for mode, smax in setup.modes_sc.items():
                if smax == 0:
                    continue
                degs = [int(x) for x in list(cst[format_name(mode)].keys())]
                if smax not in degs:
                    raise IOError
                degs = [int(d) for d in degs]
                if max(degs) < smax:
                    raise IOError

            for mode, smm in setup.modes_cc.items():
                if smm[1] == 0:
                    continue
                degs = [int(x) for x in list(cst[format_name(mode)].keys())]
                if smm[1] not in degs:
                    raise IOError
                degs = [int(d) for d in degs]
                if max(degs) < smm[1]:
                    raise IOError

            for mode, smax in setup.modes_sc_dst.items():
                if smax == 0:
                    continue
                degs = [int(x) for x in list(dst[format_name(mode)].keys())]
                if smax not in degs:
                    raise IOError
                degs = [int(d) for d in degs]
                if max(degs) < smax:
                    raise IOError

            # Is there cc for dst already implemented?
            # for mode, smm in setup.modes_cc_dst.items():
            #     if smm[1] == 0:
            #         continue
            #     degs = list(dst[format_name(mode)].keys())
            #     if str(smm[1]) not in degs:
            #         raise IOError
            #     degs = [int(d) for d in degs]
            #     if max(degs) != smm[1]:
            #         raise IOError
            return cst, dst

        except Exception as e:
            if verbose is True:
                print(e)
            else:
                pass
    if verbose is True:
        print('No database entry found: calculating')
    if bin_path is None:
        bins = ['/quanta1/home', '/net/home']
        for path in bins:
            if os.path.exists(path):
                bin_path = path
                break
        if bin_path is None:
            raise IOError('Bin_path not valid')

    cstCRUST = "%s/simons/bin/mdcplmrho_all_cstCRUST" % bin_path
    cc_cstCRUST = "%s/simons/bin/mdcplmrho_allC_cstCRUST" % bin_path
    modelfile = None

    if model == 'S20RTS':
        cstS20RTS = "%s/talavera/bin/mdcplmrho_all_cstS20RTS" % bin_path
        cc_cstS20RTS = "%s/talavera/bin/mdcplmrho_allC_cstS20RTS" % bin_path
        dstS20RTS = "%s/talavera/bin/mdcplmrho_all_dstS20RTS" % bin_path
        cc_dstS20RTS = "%s/talavera/bin/mdcplmrho_allC_dstS20RTS" % bin_path
        _maxmdeg = 20 # cst model
        _maxcdeg = 20 # crust model
        _maxddeg = 20 # dst model

    if model == 'S40RTS' or model == 'QRFSI12':
        cstS20RTS = "%s/talavera/bin/mdcplmrho_all_cstS40RTS" % bin_path
        cc_cstS20RTS = "%s/talavera/bin/mdcplmrho_allC_cstS40RTS" % bin_path
        sphm = "%s/talavera/dta/QRFSI12_l9.in" % bin_path
        dstS20RTS = "%s/talavera/bin/mdcplmrho_all_dstQRFSI12" % bin_path
        cc_dstS20RTS = "%s/talavera/bin/mdcplmrho_allC_dstQRFSI12" % bin_path
        _maxmdeg = 40 # cst model
        _maxcdeg = 20 # crust model
        _maxddeg = 12 # dst model

    if model == 'SP12RTS':
        cstS20RTS = "%s/talavera/bin/mdcplmrho_all_cstSP12RTS_E" % bin_path
        cc_cstS20RTS = "%s/talavera/bin/mdcplmrho_allC_cstSP12RTS_E" % bin_path
        dstS20RTS = "%s/talavera/bin/mdcplmrho_all_dstS20RTS" % bin_path
        cc_dstS20RTS = "%s/talavera/bin/mdcplmrho_allC_dstS20RTS" % bin_path
        _maxmdeg = 12 # cst model
        _maxcdeg = 12 # crust model
        _maxddeg = 12 # dst model

    if model == 'CRUST':
        cstS20RTS = None
        cc_cstS20RTS = None
        dstS20RTS = None
        cc_dstS20RTS = None
        _maxmdeg = 20 # cst model
        _maxcdeg = 20 # crust model
        _maxddeg = 20 # dst model

    if model in ('VSXI', 'VS'):
        cstS20RTS = "{}/simons/bin/mdcplmrho_all_cst{}".format(bin_path, model)
        cc_cstS20RTS = "{}/simons/bin/mdcplmrho_allC_cst{}".format(bin_path, model)
        dstS20RTS = None
        cc_dstS20RTS = None
        _maxmdeg = 8 # cst model
        _maxcdeg = 8 # crust model
        _maxddeg = 8 # dst model

    if mdcplbin is not None:
        cstS20RTS = mdcplbin
        cc_cstS20RTS = mdcplccbin
        dstS20RTS = None
        cc_dstS20RTS = None
        _maxmdeg = 8 # cst model
        _maxcdeg = 8 # crust model
        _maxddeg = 8 # dst model
        modelfile = model

    sc_modes, cc_modes = get_mode_names(modesin, modes_ccin)
    sc_cdeg, sc_ddeg, cc_cdeg, cc_ddeg = get_mode_deg(modesin, modes_ccin)

    """
    To do:
        Read first nmpy/data/S20RTS/S20RTS_CRUST.json
        Check if modesin, modes_ccin are alread in json
        If needed calculate missing degrees, add to json file
    """

    # Create temporary folder in which all calculations are done
    cwd = os.getcwd()
    tmp_path = 'tmp_S20'
    if os.path.exists(tmp_path):
        tmp_path = find_unique_name(tmp_path)
    os.makedirs(tmp_path)
    os.chdir(tmp_path)

    sc_coeff = {}
    if cc_modes is not None:
        if len(cc_modes) != 0:
            cc_coeff = {}

    count = 0
    for mode, s_max in zip(sc_modes, sc_cdeg):
        m = split_digit_nondigit(mode)
        if int(s_max) == 0 and int(m[2]) > 0:
            continue
        count += 1

        if (int(s_max) > 0 and int(m[2]) > 0) or (int(m[2]) == 0):
            os.system('echo "1" > modes.in')
            os.system('echo "%03d %s %03d" >> modes.in' % (int(m[0]),
                                                           m[1].lower(),
                                                           int(m[2])))
            if model not in ('SP12RTS', 'CRUST'):
                # S20RTS prediction
                for s in np.arange(0, int(s_max)+1, 2):
                    # only input coupling degrees, No degree higher then 20
                    if s not in max_sc_degrees(int(m[2])) or s > _maxmdeg:
                        continue
                    if modelfile is not None:
                        if type(modelfile) != list:
                            os.system('echo "{}" > input'.format(modelfile))
                            os.system('echo "%s" >> input' % s)
                        else:
                            os.system('echo "{}" > input'.format(modelfile[0]))
                            os.system('echo "{}" >> input'.format(modelfile[1]))
                            os.system('echo "%s" >> input' % s)
                    else:
                        os.system('echo "%s" > input' % s)
                    res = subprocess.Popen('%s < input' % cstS20RTS, shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

                    output, error = res.communicate()
                    os.system('cat mcst.dat >> mcst-S20RTS.dat')
                    os.remove('mcst.dat')
                with open('mcst-S20RTS.dat', 'r') as fh_s20rts:
                    c_s20rts_tmp = np.genfromtxt(fh_s20rts)
                if int(s_max) > _maxmdeg:
                    # Fill everything bigger than 20 with PREM values
                    for _s in range(_maxmdeg+2, int(s_max)+2, 2):
                        sdiff = (2*_s)+1
                        c_s20rts_tmp = np.hstack((c_s20rts_tmp, np.zeros(sdiff)))
                os.remove('mcst-S20RTS.dat')
            elif model == 'SP12RTS':
                # SP12RTS vs prediction
                for s in np.arange(0, int(s_max)+1, 2):
                    # only input coupling degrees, No degree higher then 20
                    if s not in max_sc_degrees(int(m[2])) or s > _maxmdeg:
                        continue
                    os.system('echo "%s" > input' % s)
                    res = subprocess.Popen('%s%s < input' % (cstS20RTS, "S"),
                                           shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

                    output, error = res.communicate()
                    os.system('cat mcst.dat >> mcst-S20RTS.dat')
                    os.remove('mcst.dat')
                with open('mcst-S20RTS.dat', 'r') as fh_s20rts:
                    c_vs_tmp = np.genfromtxt(fh_s20rts)
                    #c_s20rts_tmp = np.genfromtxt(fh_s20rts)
                os.remove('mcst-S20RTS.dat')

                # SP12RTS vp prediction
                for s in np.arange(0, int(s_max)+1, 2):
                    # only input coupling degrees, No degree higher then 20
                    if s not in max_sc_degrees(int(m[2])) or s > _maxmdeg:
                        continue
                    os.system('echo "%s" > input' % s)
                    res = subprocess.Popen('%s%s < input' % (cstS20RTS, "P"),
                                           shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

                    output, error = res.communicate()
                    os.system('cat mcst.dat >> mcst-S20RTS.dat')
                    os.remove('mcst.dat')
                with open('mcst-S20RTS.dat', 'r') as fh_s20rts:
                    c_vp_tmp = np.genfromtxt(fh_s20rts)
                    #c_s20rts_tmp = np.genfromtxt(fh_s20rts)
                os.remove('mcst-S20RTS.dat')

                if int(s_max) > _maxmdeg:
                    # Fill everything bigger than 20 with PREM values
                    for _s in range(_maxmdeg+2, int(s_max)+2, 2):
                        sdiff = (2*_s)+1
                        c_vs_tmp = np.hstack((c_vs_tmp, np.zeros(sdiff)))
                        c_vp_tmp = np.hstack((c_vp_tmp, np.zeros(sdiff)))
                        #c_s20rts_tmp = np.hstack((c_s20rts_tmp, np.zeros(sdiff)))

                c_s20rts_tmp = np.add(c_vs_tmp, c_vp_tmp)

            # CRUST prediction
            if include_CRUST is True:
                for s in np.arange(0, int(s_max)+1, 2):

                    # only input coupling degrees, No degree higher then 20
                    if s not in max_sc_degrees(int(m[2])) or s > _maxcdeg:
                        continue
                    os.system('echo "%s" > input' % s)
                    res = subprocess.Popen('%s < input' % cstCRUST, shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

                    output, error = res.communicate()
                    os.system('cat mcst.dat >> mcst-CRUST.dat')
                    os.remove('mcst.dat')
                with open('mcst-CRUST.dat', 'r') as fh_crust:
                    c_crust_tmp = np.genfromtxt(fh_crust)
                if int(s_max) > _maxcdeg:
                    # Fill everything bigger than 20 with PREM values
                    for _s in range(_maxcdeg+2, int(s_max)+2, 2):
                        sdiff = (2*_s)+1
                        c_crust_tmp = np.hstack((c_crust_tmp, np.zeros(sdiff)))

                os.remove('mcst-CRUST.dat')
                if model == 'CRUST':
                    sc_coeff[mode] = c_crust_tmp.transpose()
                else:
                    sc_coeff[mode] = np.add(c_s20rts_tmp, c_crust_tmp).transpose()
            else:
                sc_coeff[mode] = c_s20rts_tmp.transpose()
            # Cross coupling coefficients
            # CC after each mode in Self
            # 1-1  ->  1-2  ->  1-3
            #          2-2  ->  2-3
            #                   3-3
            if cc_modes is not None:
                if len(cc_modes) == 0:
                    continue
                # Workaround for cross-coupling
                if setup is not None:
                    cc_cdeg = list(setup.modes_cc.values())
                else:
                    # ??????
                    x = [list(setup.modes_cc.values())[0][0]]
                    cc_cdeg = x + [int(cc_cdeg[0])]
                    cc_cdeg = [cc_cdeg]

                # reading for 2 and 3 modes
                if count == 1 and len(cc_modes) >= 1:
                    # 1-2  ->  1-3
                    cc_modesin = cc_modes[0:2]
                    cc_cdegin = cc_cdeg[0:2]

                elif count == 2 and len(cc_modes) == 3:
                    # ->  2-3
                    cc_modesin = cc_modes[-1:]
                    cc_cdegin = cc_cdeg[-1:]
                else:
                    cc_modesin = []
                    cc_cdegin = []

                for modecc, sdeg in zip(cc_modesin, cc_cdegin):
                    if sdeg[1] > 0:
                        cc = modecc.split('-')
                        max_cc_input = []  # only calculate cpl s degrees

                        # number of cpl modes is always 2
                        os.system('echo "%s" > modes.in' % 2)
                        for i, m in enumerate(cc):
                            m = split_digit_nondigit(m)
                            os.system('echo "%03d %s %03d" >> modes.in'
                                      % (int(m[0]), m[1].lower(), int(m[2])))
                            for _m in m:
                                max_cc_input += [_m]

                        ccdegs = max_cc_degrees(max_cc_input)

                        # S20RTS cc prediction
                        if model != 'CRUST':
                            for s in np.arange(sdeg[0], int(sdeg[1])+1, 2):
                                # only input coupling degrees
                                if s not in ccdegs:
                                    continue
                                if os.path.exists('input'):
                                    os.remove('input')
                                if modelfile is not None:
                                    if type(modelfile) != list:
                                        os.system('echo "{}" > input'.format(modelfile))
                                        os.system('echo "%s" >> input' % s)
                                    else:
                                        os.system('echo "{}" > input'.format(modelfile[0]))
                                        os.system('echo "{}" >> input'.format(modelfile[1]))
                                        os.system('echo "%s" >> input' % s)
                                else:
                                    os.system('echo "%s" > input' % s)
                                res = subprocess.Popen('%s < input' % cc_cstS20RTS,
                                                       shell=True,
                                                       stdout=subprocess.PIPE,
                                                       stderr=subprocess.PIPE)

                                output, error = res.communicate()
                                os.system('cat mcst.dat >> mcst-S20RTS.dat')
                                os.remove('mcst.dat')
                            with open('mcst-S20RTS.dat', 'r') as fh:
                                c_s20rts_tmp = np.genfromtxt(fh)
                            os.remove('mcst-S20RTS.dat')
                        # CRUST cc prediction
                        if include_CRUST is True:
                            for s in np.arange(sdeg[0], int(sdeg[1])+1, 2):
                                # only input coupling degrees
                                if s not in ccdegs:
                                    continue
                                os.system('echo "%s" > input' % s)
                                res = subprocess.Popen('%s < input' % cc_cstCRUST,
                                                       shell=True,
                                                       stdout=subprocess.PIPE,
                                                       stderr=subprocess.PIPE)

                                output, error = res.communicate()
                                os.system('cat mcst.dat >> mcst-CRUST.dat')
                                os.remove('mcst.dat')

                            with open('mcst-CRUST.dat', 'r') as fh:
                                c_crust_tmp = np.genfromtxt(fh)
                            os.remove('mcst-CRUST.dat')
                            if model == 'CRUST':
                                cc_coeff[modecc] = c_crust_tmp.transpose()
                            else:
                                cc_coeff[modecc] = np.add(c_s20rts_tmp,
                                                          c_crust_tmp).transpose()
                        else:
                            cc_coeff[modecc] = c_s20rts_tmp.transpose()
                        # os.system('cat mcst.dat >> cst.dat')
                        # os.remove('mcst.dat')
            else:
                cc_modesin = None

    # dst calculation here?
    if modes_dst is None:
        dst = None
        sc_coeff_dst = None
    else:
        sc_ddeg = get_mode_deg_dst(modes_dst)
        sc_modes_dst, cc_modes_dst = get_mode_names(modes_dst, modes_cc_dst)
        sc_coeff_dst = {}
        cc_coeff_dst = {}
        count = 0
        for mode, s_max in zip(sc_modes_dst, sc_ddeg):
            m = split_digit_nondigit(mode)
            count += 1
            if int(s_max) >= 2:
                os.system('echo "1" > modes.in')
                os.system('echo "%03d %s %03d" >> modes.in' % (int(m[0]),
                                                               m[1].lower(),
                                                               int(m[2])))
                for s in np.arange(2, int(s_max)+1, 2):
                    # only input coupling degrees
                    if s not in max_sc_degrees(int(m[2]))[1:] or s > _maxddeg:
                        continue
                    # S20RTS prediction
                    if model == "S20RTS" or model == 'SP12RTS':
                        os.system('echo %s > input' % R)  # R=-0.2 wrt vs
                        os.system('echo "%s" >> input' % s)
                        res = subprocess.Popen('%s < input' % dstS20RTS,
                                               shell=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
                    elif model == 'S40RTS' or model == 'QRFSI12':
                        os.system('echo %s > input' % s)
                        res = subprocess.Popen('%s -m %s < input' % (dstS20RTS, sphm),
                                               shell=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
                    output, error = res.communicate()
                os.system('cat dst.dat >> mcst-S20RTS.dat')
                os.remove('dst.dat')
                with open('mcst-S20RTS.dat', 'r') as fh_s20rts:
                    c_s20rts_tmp = np.genfromtxt(fh_s20rts)
                    os.remove('mcst-S20RTS.dat')

                if int(s_max) > _maxddeg:
                    # Fill everything bigger than model max with PREM values
                    for _s in range(_maxddeg+2, int(s_max)+2, 2):
                        sdiff = (2*_s)+1
                        c_s20rts_tmp = np.hstack((c_s20rts_tmp, np.zeros(sdiff)))

                sc_coeff_dst[mode] = c_s20rts_tmp.transpose()

                # Cross coupling coefficients
                # CC after each mode in Self
                # 1-1  ->  1-2  ->  1-3
                #          2-2  ->  2-3
                #                   3-3
                if modes_cc_dst is not None:
                    # Workaround for cross-coupling
                    if setup is not None:
                        cc_cdeg = list(setup.modes_cc_dst.values())
                    else:
                        # ??????
                        x = [list(setup.modes_cc_dst.values())[0][0]]
                        cc_cdeg = x + [int(cc_cdeg[0])]
                        cc_cdeg = [cc_cdeg]

                    # reading for 2 and 3 modes
                    if count == 1 and len(cc_modes_dst) >= 1:
                        # 1-2  ->  1-3
                        cc_modesin = cc_modes_dst[0:2]
                        cc_cdegin = cc_cdeg[0:2]

                    elif count == 2 and len(cc_modes_dst) == 3:
                        # ->  2-3
                        cc_modesin = cc_modes_dst[-1:]
                        cc_cdegin = cc_cdeg[-1:]
                    else:
                        cc_modesin = []
                        cc_cdegin = []

                    for modecc, sdeg in zip(cc_modesin, cc_cdegin):
                        if sdeg[1] > 0:
                            cc = modecc.split('-')
                            max_cc_input = []  # only calculate cpl s degrees

                            # number of cpl modes is always 2
                            os.system('echo "%s" > modes.in' % 2)
                            for i, m in enumerate(cc):
                                m = split_digit_nondigit(m)
                                os.system('echo "%03d %s %03d" >> modes.in'
                                          % (int(m[0]), m[1].lower(), int(m[2])))
                                for _m in m:
                                    max_cc_input += [_m]

                            ccdegs = max_cc_degrees(max_cc_input)

                            # S20RTS cc prediction
                            for s in np.arange(sdeg[0], int(sdeg[1])+1, 2):
                                # only input coupling degrees
                                if s not in ccdegs:
                                    continue

                                # S20RTS prediction
                                if model == "S20RTS" or model == 'SP12RTS':
                                    os.system('echo %s > input' % R)  # R=-0.2 wrt vs
                                    os.system('echo "%s" >> input' % s)
                                    res = subprocess.Popen('%s < input' % cc_dstS20RTS,
                                                           shell=True,
                                                           stdout=subprocess.PIPE,
                                                           stderr=subprocess.PIPE)
                                elif model == 'S40RTS' or model == 'QRFSI12':
                                    os.system('echo %s > input' % s)
                                    res = subprocess.Popen('%s -m %s < input' % (cc_dstS20RTS, sphm),
                                                           shell=True,
                                                           stdout=subprocess.PIPE,
                                                           stderr=subprocess.PIPE)

                                output, error = res.communicate()
                                os.system('cat dst.dat >> mcst-S20RTS.dat')
                                os.remove('dst.dat')
                            with open('mcst-S20RTS.dat', 'r') as fh:
                                c_s20rts_tmp = np.genfromtxt(fh)
                            os.remove('mcst-S20RTS.dat')
                            cc_coeff_dst[modecc] = c_s20rts_tmp.transpose()
                            # os.system('cat mcst.dat >> cst.dat')
                            # os.remove('mcst.dat')
                else:
                    cc_coeff_dst = None

    for _file in ('modes.in', 'cst.dat', 'mdcpl.out', 'raw.dat', 'input'):
        try:
            # os.system('cp mcst.dat cst.dat')
            os.remove(_file)
        except FileNotFoundError:
            if verbose is True:
                print(_file)
            pass

    cT = sc_coeff[sc_modes[0]]
    mode = sc_modes[0]
    count = 1
    for mode in sc_modes[1:]:
        if cc_modes is not None:
            if len(cc_modes) == 0:
                continue
            # reading for 2 and 3 modes
            if count == 1 and len(cc_modes) >= 1:
                # 1-2  ->  1-3
                cc_modesin = cc_modes[0:2]
                cc_cdegin = cc_cdeg[0:2]
            elif count == 2 and len(cc_modes) == 3:
                # ->  2-3
                cc_modesin = cc_modes[-1:]
                cc_cdegin = [cc_cdeg[-1]]
            else:
                cc_modesin = []
                cc_cdegin = []

            for modecc, sdeg in zip(cc_modesin, cc_cdegin):
                if sdeg[1] > 0:
                    cT = np.hstack((cT, cc_coeff[modecc]))
        if mode in sc_coeff:
            cT = np.hstack((cT, sc_coeff[mode]))
        count += 1

    # Only works for 2 modes, 1 measured and 1 as PREM
    # when the CC between them is also measured
    if not sc_modes[1:] and cc_modesin:
        for modecc, sdeg in zip(cc_modesin, cc_cdegin):
            if sdeg[1] > 0:
                cT = np.hstack((cT, cc_coeff[modecc]))

    count = 0
    if sc_coeff_dst:
        for mode, s_max in zip(sc_modes_dst, sc_ddeg):
            count += 1
            cT = np.hstack((cT, sc_coeff_dst[mode]))

            if modes_cc_dst is not None:
                # reading for 2 and 3 modes
                if count == 1 and len(cc_modes_dst) >= 1:
                    # 1-2  ->  1-3
                    cc_modesin = cc_modes_dst[0:2]
                    cc_cdegin = cc_cdeg[0:2]
                elif count == 2 and len(cc_modes_dst) == 3:
                    # ->  2-3
                    cc_modesin = cc_modes_dst[-1:]
                    cc_cdegin = [cc_cdeg[-1]]
                else:
                    cc_modesin = []
                    cc_cdegin = []
                for modecc, sdeg in zip(cc_modesin, cc_cdegin):
                    if sdeg[1] > 0:
                        cT = np.hstack((cT, cc_coeff_dst[modecc]))
    # Preparing cst and dst files
    if keep_mcst:
        # d00 needs to be converted to haydars format
        with open('mcst_zero.dat', 'w') as f:
            for c in cT:
                f.write("%s\n" % c)
    cst, dst = get_cst(modes=modesin, modes_cc=modes_ccin, c=cT, noc=False,
                       modes_dst=modes_dst, modes_cc_dst=modes_cc_dst,)
    #print(modes_ccin)
    # scfiles = glob.glob("mcst-*_sc_*.dat")
    # ccfiles = glob.glob("mcst-*_cc_*.dat")
    # for f in scfiles:
    #     os.remove(f)
    # for f in ccfiles:
    #     os.remove(f)

    # Switching back to work directory and remove tmp directory
    os.chdir(cwd)
    os.removedirs(tmp_path)

    # Only write the default setting to db
    # if R == -0.2 and model != 'QRFSI12':
    WRITE2DB = False

    if model in ('S20RTS', 'S40RTS'):
        WRITE2DB = True
    elif R == -0.2:
        WRITE2DB = True
    if mdcplbin is not None:
        WRITE2DB = False

    if WRITE2DB is True:
        _write_cst_S20RTS_db(cst, dst, file_name, verbose=verbose)
    return cst, dst


def _read_pickle(filename, **kwargs):
    """
    Read and return Modes from pickled Modes file.

    .. warning::
        This function should NOT be called directly, it registers via the
        nmPy :func:`~frospy.core.modes.read` function, call this instead.

    :type filename: str
    :param filename: Name of the pickled Modes file to be read.
    :rtype: :class:`~frospy.core.modes.Modes`
    :return: A Modes object.
    """
    kwargs = {}

    if isinstance(filename, (str, native_str)):
        with open(filename, 'rb') as fp:
            return pickle.load(fp, **kwargs)
    else:
        return pickle.load(filename, **kwargs)


def get_modes4cst(modes):
    if type(modes) != list:
        modes = [modes]
    modes = [x.upper() for x in modes]
    modesin = ['-']
    modes_ccin = ['-']
    for m in modes:
        if '-' in m:
            _cc = ''
            for _m in m.split('-'):
                _cc += ' ' + ' '.join(split_digit_nondigit(_m))
                _sc = "%s 20 0" % ' '.join(split_digit_nondigit(_m))
                if _sc not in modesin:
                    modesin += [_sc]
            modes_ccin += [_cc + ' - -']

        else:
            _sc = ' '.join(split_digit_nondigit(m)) + ' 20 0'
            if _sc not in modesin:
                modesin += [_sc]
    modes_ccin = [x.upper() for x in modes_ccin]
    modes_ccin[0] = str(len(modes_ccin) - 1)
    modesin = [x.upper() for x in modesin]
    modesin[0] = str(len(modesin) - 1)

    if modes_ccin == ['-']:
        modes_ccin = None
    if modesin == ['-']:
        modesin = None

    sc, cc = get_mode_names(modesin, modes_ccin)
    modes_sc = Modes()
    modes_cc = Modes()
    modes_all = read_modes()
    for m in sc:
        modes_sc += modes_all.select(name=m)

    if len(modes) == 1:
        cc = None
        modes_cc = None
        modes_ccin = None

    if cc is not None:
        for m in cc:
            header = {'n': -1, 'type': 'CC', 'l': -1, 'name': m,
                      'sens': None, 'freq': 0, 'Q': 0}
            modes_cc += Mode(header)

    return modes_sc, modes_cc, modesin, modes_ccin
