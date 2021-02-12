from __future__ import absolute_import, print_function

from obspy.core.util import AttribDict
from frospy.core.modes import format_name
from frospy.core.splittingfunc.splittingfunc import get_model_cst

from frospy.splitting.load import loadmodel
from frospy.util.base import N_splitting_coeffs, split_digit_nondigit

import numpy as np

from collections import OrderedDict


def buildmodel(setup):
    if (
        setup.startmodel in ('PREM', 'S20RTS', 'AD', 'RR') or
        type(setup.startmodel) in (dict, OrderedDict, AttribDict)
       ):
        # Check if modes are actually measured
        if type(setup.startmodel) in (dict, OrderedDict, AttribDict):
            ModeErr = False
            msg = ""
            for mode in setup.startmodel.keys():
                if mode not in setup.modes_sc.keys():
                    if mode not in setup.modes_cc.keys():
                        _msg = 'startmodel mode: {} not in setup\n'
                        msg += _msg.format(mode)
                        ModeErr = True
            if ModeErr is True:
                raise IOError(msg)

        # maxmdeg sets the maximum degree for model
        # setup.startmodel should be list or dict-like, with the name of the
        # model as first entry / key and the maximum model degree as a second
        # entry / value. The remaining degrees will be filled with zeros
        model_modes = None
        if type(setup.startmodel) in (dict, OrderedDict, AttribDict):
            if sorted(list(setup.startmodel.keys())) == ['cst', 'dst']:
                mzero = OrderedDict([('cst', setup.startmodel['cst']),
                                     ('dst', setup.startmodel['dst'])])
                return mzero
            maxmdeg = list(setup.startmodel.keys())
            model = list(setup.startmodel.values())

            model_modes = [format_name(x) for x in maxmdeg]
            maxmdeg = None
        else:
            maxmdeg = None
            model = setup.startmodel
        mzero = OrderedDict()
        mzero_cst = []
        mzero_dst = None
        countsc = 1
        countcc = 1
        for mode, smax in setup.modes_sc.items():
            ldeg = int(split_digit_nondigit(mode)[-1])
            # radial modes and other modes with smax > 0
            if model_modes is not None:
                if format_name(mode) in model_modes:
                    model = setup.startmodel[mode]
                else:
                    model = 'PREM'
            if ((smax > 0 and ldeg > 0) or (ldeg == 0)):
                coeffs = get_model_cst(setup, model, mode, 0, smax)
                mzero_cst.append((mode, coeffs))

            if len(setup.modes_cc) > 0:
                for modecc, s in setup.modes_cc.items():
                    if model_modes is not None:
                        if format_name(modecc) in model_modes:
                            model = setup.startmodel[mode]
                        else:
                            model = 'PREM'
                    smin = s[0]
                    if maxmdeg is not None:
                        smax = maxmdeg
                    else:
                        smax = s[1]

                    # CC with smax > 0
                    if smax > 0:
                        if countsc == 1 and countcc > 2:
                            pass
                        elif countsc == 3:
                            pass
                        coeffs = get_model_cst(setup, model, modecc,
                                               smin, smax)
                        if maxmdeg is not None:
                            Z = np.zeros(N_splitting_coeffs(maxmdeg, s[1]))
                            coeffs = np.hstack((coeffs, Z))
                        else:
                            mzero_cst.append((modecc, coeffs))
                    countcc += 1
            countsc += 1
            countcc = 1
        cst = OrderedDict(mzero_cst)
        mzero.update(OrderedDict([("cst", cst)]))
        if len(setup.modes_sc_dst) > 0:
            mzero_dst = []
            for mode, smax in setup.modes_sc_dst.items():
                if model_modes is not None:
                    if format_name(mode) in model_modes:
                        model = setup.startmodel[mode]
                    else:
                        model = 'PREM'

                # in dst: smax > 2
                if smax >= 2:
                    if (format_name(mode) in model_modes):
                        model = setup.startmodel[mode]
                    else:
                        model = 'PREM'
                    coeffs = get_model_cst(setup, model, mode, 2,
                                           smax, kind="dst")
                    mzero_dst.append((mode, coeffs))

            dst = OrderedDict(mzero_dst)
            mzero.update(OrderedDict([("dst", dst)]))

    elif setup.startmodel is None:
        mzero = OrderedDict([('cst', None), ('dst', None)])
    elif type(setup.startmodel) in (np.ndarray, list):
        mzero = OrderedDict([('cst', setup.startmodel[0]),
                             ('dst', setup.startmodel[1])])
    else:
        try:
            # if isfile(setup.startmodel):
            # Needs checking that the model is the correct size
            print('Loading mcst model from file')
            mzero = loadmodel(setup=setup, ifile=setup.startmodel,
                              name='', damp='')
        except Exception:
            print('Invalid input, no mzero.dat file will be created')
            mzero = None
    return mzero


def build_deg_dependent_damp(setup, modes='all', degs='all', function='linear',
                             smag_diff=10):
    def get_m(x_min, x_max, function, x, y_max=smag_diff-1):
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
            a = y_max / np.log(x_max+1)
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
