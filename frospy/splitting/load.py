from __future__ import absolute_import, print_function

from frospy.core.splittingfunc.read import (read_cst, _read_pickle,
                                            get_modes4cst)
from frospy.core.splittingfunc.splittingfunc import SplittingFunc
from frospy.core.splittingfunc.splittingfunc import get_header
from frospy.core.splittingfunc.set import Set
from frospy.core.modes import format_name
from frospy.core.database.query import db_query
from frospy.core.setup.settings import read as read_setup
import os


def loadmodel(*args):
    return load(*args)


def load(ifile=None, modes=None, setup=None, modesin_dir=None,
         format=None, name='data', damp=None, R=-0.2, db_model=None,
         return_set=False, include_CRUST=True,
         verbose=False, name_overide=False,
         mdcplbin=None, mdcplccbin=None):

    """
    param setup: :frospy.core.setup.settings.Setup object:
    param ifile: path to file
    param modesin_dir: path to folder containing modes.in and modes_cc.in
    param format: 'S20RTS', 'REM', 'RR', 'TZ' or 'pickle'
    """

    models = ['S20RTS', 'S40RTS', 'REM', 'RR', 'TZ', 'CB', 'TCB', 'AD', 'PREM',
              'HT', 'QM1', 'DE', 'GLW', 'GD', 'PK', 'MW', 'WZM', 'SAS',
              'STS_SC', 'STS_GC_SC', 'STS_GC_CC', 'Sumatra',
              'S20RTS+CRUST+BT', 'S20RTS+CRUST+Tr',
              'S20RTS+CRUST+Wh', 'S20RTS+CRUST+Ro',
              'SP12RTS', 'QRFSI12', 'CRUST', 'VSXI', 'VS',
              'BT', 'Tr', 'Ro', 'Wh',
              '$\phi$=1.20', '$\phi$=1.10', '$\phi$=1.04',
              '$\phi$=0.96', '$\phi$=0.90', '$\phi$=0.80',
              'custom']

    model = None
    if modes is not None:
        if type(modes) == str:
            modes = [modes]
        for i, _m in enumerate(modes):
            modes[i] = format_name(_m)

    # WE HAVE TO COMMENT ON THE IF CONDITIONS HERE!!!

    # if setup is None and ifile is not None:
    #     if not ifile.endswith('.sqlite3'):
    #         msg = "Error in setup: Give setup object or path to setup.pickle"
    #         raise IOError(msg)

    if ifile is None and format is None:
        msg = "Error in ifile: Give path to cst-file (mcst.dat)\n"
        msg += "or\nGive format ('S20RTS', 'REM', 'RR', 'TZ')"
        raise IOError(msg)

    # if damp is None and format not in models:
    #     msg = "Error in damp: Give damping value"
    #     raise IOError(msg)

    if type(setup) == str:
        setup = read_setup(setup)
        # setup = read_setup(setup)

    if format is None and ifile is not None and not type(ifile) == list:
        if not ifile.endswith('sph'):
            # try to guess format from file extension
            _, format = os.path.splitext(ifile)
            format = format[1:]

    if format == 'pickle':
        cst, dst = _read_pickle(ifile)
        pass

    if ifile is not None and ifile.endswith('.sqlite3') and not type(ifile) == list:
        if not ifile.endswith('sph'):
            if name_overide is True:
                name = name
            else:
                name = db_model

            cst_out = read_cst(setup=setup, modes=modes, cfile=ifile,
                               model=db_model)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            if verbose is True:
                print('cst errors: \n', cst_errors)
                print('cst errors: \n', dst_errors)
                print('\n')

            if setup is not None:
                header = get_header(setup.rundir, modes_sc, modes_cc,
                                    name=name, damp=None)
            else:
                if type(modes) is not list:
                    modes = [modes]
                for _i, _name in enumerate(modes):
                    modes[_i] = format_name(_name).upper()

                damp = db_query(db_path=ifile, table=db_model, select='damp',
                                condition_dict={'mode': modes[0], 'kind': 'cst'})
                if verbose is True:
                    print('damping found: ', damp)
                header = get_header(None, modes_sc, modes_cc,
                                    name=name, damp=damp[0][0])

    elif setup is not None and ifile is not None and not type(ifile) == list:
        if not ifile.endswith('sph'):
            cst_out = read_cst(setup=setup, cfile=ifile)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            header = get_header(setup.rundir, modes_sc, modes_cc, name=name,
                                damp=damp)

    elif format in models or mdcplbin is not None:
        if setup is not None:
            cst_out = read_cst(setup=setup, cfile=format, R=R,
                               include_CRUST=include_CRUST,
                               mdcplbin=mdcplbin)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            header = get_header(setup.rundir, modes_sc, modes_cc, damp=0,
                                name=name, model=format)
        elif modesin_dir is not None:
            cst_out = read_cst(format, modesin_dir,
                               include_CRUST=include_CRUST,
                               mdcplbin=mdcplbin)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            name = format
            model = format
            header = get_header(modesin_dir, modes_sc, modes_cc,
                                name=name, model=model, damp=damp)
        elif mdcplbin is not None:
            cst_out = read_cst(cfile=ifile, modes=modes, verbose=verbose,
                               include_CRUST=include_CRUST,
                               mdcplbin=mdcplbin,
                               mdcplccbin=mdcplccbin)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            name = format
            model = format
            header = get_header(None, modes_sc, modes_cc,
                                name=name, model=model, damp=damp)
        else:
            cst_out = read_cst(cfile=format, modes=modes, verbose=verbose,
                               include_CRUST=include_CRUST)
            cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
            name = format
            model = format
            header = get_header(None, modes_sc, modes_cc,
                                name=name, model=model, damp=damp)
    elif format == 'dat':
        if modesin_dir is not None:
            cst_out = read_cst(cfile=ifile, modes_dir=modesin_dir)

        else:
            cst_out = read_cst(cfile=ifile, modes=modes)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        header = get_header(modesin_dir, modes_sc, modes_cc,
                            name=name, model=model, damp=damp)

    # if mode is not defined, it will load all modes from the file,
    # in this case we need to loop over them
    if return_set is True:
        S = Set()
        for _m, _cst in cst.items():
            c = {_m: _cst}
            d = {_m: dst[_m]}

            try:
                c_err = {_m: cst_errors[_m]}
                d_err = {_m: dst_errors[_m]}
            except Exception:
                pass
            modes_sc, modes_cc, modesin, modes_ccin = get_modes4cst(_m)

            header = get_header(modesin_dir, modes_sc,
                                modes_cc, name=name, model=model, damp=damp)
            S += SplittingFunc(header=header, cst=c, dst=d,
                               cst_errors=c_err, dst_errors=d_err)
        return S

    elif modes is not None or return_set is False:
        return SplittingFunc(header=header, cst=cst, dst=dst,
                             cst_errors=cst_errors, dst_errors=dst_errors)
