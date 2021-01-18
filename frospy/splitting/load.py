from __future__ import absolute_import, print_function

from frospy.core.splittingfunc.read import read_cst, _read_pickle
from frospy.core.splittingfunc.splittingfunc import SplittingFunc
from frospy.core.splittingfunc.splittingfunc import get_header
from frospy.core.modes import format_name
from frospy.core.database.query import db_query
import os

def update_SF_from_db(dbpath ='/tmp/eejit_simons/splitting/modes/cst_final.sqlite3'):
    # dbpath = '/quanta1/home/simons/splitting/modes/cst_final.sqlite3'
    paths = get_paths('all', results='paper2', host='net', checked=True, dictionary=True)
    SF = Set()


    for mode, values in paths.items():
        damp = float(values[0].split('/')[-1].split('d')[1][:-1])
        try:
            SF += loadmodel(ifile=dbpath, modes=format_name(mode).upper(),
                            name='data', damp='0', db_model='data')
        except IndexError:
            continue
    return SF


def loadmodel(modes=None, setup=None, ifile=None, modesin_dir=None,
              format=None, name=None, damp=None, R=-0.2, db_model=None,
              verbose=False):
    """
    param setup: :frospy.core.setup.settings.Setup object:
    param ifile: path to file
    param modesin_dir: path to folder containing modes.in and modes_cc.in
    param format: 'S20RTS', 'REM', 'RR', 'TZ' or 'pickle'
    """

    models = ['S20RTS', 'S40RTS', 'REM', 'RR', 'TZ', 'CB', 'TCB', 'AD', 'PREM',
              'HT', 'QM1', 'DE', 'GLW', 'GD', 'PK', 'MW', 'WZM', 'SAS', 'Sumatra',
              'S20RTS+CRUST+BT', 'S20RTS+CRUST+Tr',
              'S20RTS+CRUST+Wh', 'S20RTS+CRUST+Ro',
              'BT', 'Tr', 'Ro', 'Wh',
              '$\phi$=1.20', '$\phi$=1.10', '$\phi$=1.04',
              '$\phi$=0.96', '$\phi$=0.90', '$\phi$=0.80']

    model = None

    if modes is not None:
        if type(modes) == str:
            modes = [modes]
        for i, _m in enumerate(modes):
            modes[i] = format_name(_m)

    # WE HAVE TO COMMENT ON THE IF CONDITIONS HERE!!!

    if setup is None and ifile is not None:
        if not ifile.endswith('.sqlite3'):
            msg = "Error in setup: Give setup object or path to setup.pickle"
            raise IOError(msg)

    if ifile is None and format is None:
        msg = "Error in ifile: Give path to cst-file (mcst.dat)\n"
        msg += "or\nGive format ('S20RTS', 'REM', 'RR', 'TZ')"
        raise IOError(msg)

    # if damp is None and format not in models:
    #     msg = "Error in damp: Give damping value"
    #     raise IOError(msg)

    if type(setup) == str:
        raise IOError('setup, must be object, not file')
        # setup = read_setup(setup)

    if format is None and ifile is not None:
        # try to guess format from file extension
        _, format = os.path.splitext(ifile)
        format = format[1:]

    if format == 'pickle':
        cst, dst = _read_pickle(ifile)
        pass

    if ifile is not None and ifile.endswith('.sqlite3'):
        cst_out = read_cst(setup=setup, modes=modes, cfile=ifile,
                           model=db_model)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        if verbose is True:
            print('cst errors: \n', cst_errors)
            print('cst errors: \n', dst_errors)
            print('\n')

        if setup is not None:
            header = get_header(setup.rundir, modes_sc, modes_cc,
                                name=db_model, damp=None)
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
                                name=db_model, damp=damp[0][0])

    if setup is not None and ifile is not None:
        cst_out = read_cst(setup=setup, cfile=ifile)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        header = get_header(setup.rundir, modes_sc, modes_cc, name=name,
                            damp=damp)

    if format in models and setup is not None:
        cst_out = read_cst(setup=setup, cfile=format, R=R)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        header = get_header(setup.rundir, modes_sc, modes_cc, damp=0,
                            name=name, model=format)

    elif format == 'dat' and modesin_dir is not None:
        cst_out = read_cst(ifile, modesin_dir)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]

    elif (format in models and modesin_dir is not None):
        cst_out = read_cst(format, modesin_dir)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        name = format
        model = format

    elif format in models and modes is not None:
        cst_out = read_cst(cfile=format, modes=modes, verbose=verbose)
        cst, dst, cst_errors, dst_errors, modes_sc, modes_cc = cst_out[:]
        name = format
        model = format
        header = get_header(None, modes_sc, modes_cc,
                            name=name, model=model, damp=damp)

    else:
        header = get_header(modesin_dir, modes_sc, modes_cc,
                            name=name, model=model, damp=damp)

    splitf = SplittingFunc(header=header, cst=cst, dst=dst,
                           cst_errors=cst_errors, dst_errors=dst_errors)

    return splitf
