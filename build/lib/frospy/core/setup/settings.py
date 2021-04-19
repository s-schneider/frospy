import inspect
import os
import glob
from os.path import basename

from obspy.core import AttribDict
import numpy as np
from collections import OrderedDict
from abc import ABCMeta, abstractmethod
import json

from frospy.core.modes import read as read_modes

from frospy.util.base import (N_splitting_coeffs, split_digit_nondigit,
                              max_cc_degrees)
from frospy.util.write import write_pickle
from frospy.util.read import read_pickle


class DefaultSetup(object):
    __metaclass__ = ABCMeta
    """
    Default values for settings
    """

    # binary file locations
    bins = AttribDict()
    bins['synseis'] = 'synseis-new-ah'
    bins['synseis_inv'] = 'synseis-inv'
    bins['rot_ellip'] = 'mdcplmrho_all_rotellip'
    bins['write_dst'] = 'write-dst'
    bins['mdcpl'] = 'matrixcst'
    bins['mdcplmq'] = 'matrixdst'
    bins['addmdcplmq'] = 'addmdcplmq'  # adds matrixcst and matrixdst
    bins['matdiag'] = 'matdiag'
    bins['cst_partials'] = 'compu-deriv-cst'
    bins['dst_partials'] = 'compu-deriv-dst'
    bins['fsi_partials'] = 'compu-deriv-inv'
    bins['buildATA'] = 'buildATA'
    bins['addATA'] = 'merge.py'
    bins['invATA'] = 'invATA-new'
    bins['avmisfit'] = 'avmisfit'
    bins['avmisfit_allmodes'] = 'avmisfit-allmodes'
    bins['posterioriCd'] = 'posterioriCd'

    bindir = '/home/deuss/bin'
    rundir = None
    inversion_outdir = 'inversion_out'
    default_data = {
        'stig': '//nfs/stig/deuss/modes/alldatafiles',
        'quanta': '/quanta1/home/simons/splitting/alldatafiles',
        'scratch': '/scratch/simons/splitting/alldatafiles',
        'local': '/home/simons/dev/alldatafiles'
        }
    nmpydir = ''

    channel = 'VH'

    use_Qfiles = False

    # Specify pbs
    scripttype = 'slurm'
    pbsjob = 'test'
    walltime = None
    partition = 'allq'
    nofcpu = 48
    nodes = 1
    exclude = None

    iterno = 1
    damping = 0.01
    startmodel = None  # 'PREM'
    dst_R = -0.2
    rotation = 1
    ellip = 1
    modes_cc = OrderedDict()
    modes_sc_dst = OrderedDict()
    segmentsdir = AttribDict()
    modes = None
    comp = None
    keep_all_output = True
    seg_suffix = 'dat'
    schema = None
    uncertainties = AttribDict({'files': None, 'damping': None})
    processing = []

    @abstractmethod
    def __init__(self, args):
        pass


class Setups(AttribDict):
    def __init__(self, intype, args=None):
        pass

    def __str__(self, extended=False):
        _pstr = "Setup for "
        if self.intype is not None:
            _pstr += "%s " % self.intype
        _pstr += "Inversion \n\n"
        return _pstr

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))


class Setup(DefaultSetup):

    def __init__(self, intype, args=None, sanity_checks=True):
        if intype == 'CST':
            self.intype = intype
            self.run_keys = ['mdcpl', 'matdiag', 'synseis', 'cstpar',
                             'cstsolver']
        if intype == 'FSI':
            self.intype = intype
            self.run_keys = ['mdcpl', 'matdiag', 'synseis', 'fsipar',
                             'fsisolver']
        if intype == 'synthetics':
            self.intype = intype
            self.run_keys = ['mdcpl', 'matdiag', 'synseis']

        # Set defaults
        if args is not None:
            for key, value in args.items():
                if key == 'modes' and intype == 'CST':
                    modes_sanity_check(key, value)
                    # Create list of names
                    m_list = get_mode_list(value)
                    self.modes = read_modes(modenames=m_list)
                    self.modes_sc = value
                elif key == 'modes_sc_dst' and intype == 'CST':
                    if 'modes' not in args:
                        raise IOError("You need modes to use modes_dst")
                        if sanity_checks is True:
                            modes_sanity_check(key, value)
                    self.modes_sc_dst = value
                elif key == 'modes_cc' and intype == 'CST':
                    modes_sanity_check(key, value)
                    self.modes_cc = value
                    self.bins['mdcpl'] = 'matrixcstC'
                    self.bins['cst_partials'] = 'compu-deriv-cstC'
                elif key == 'damping':
                    self.damping = get_value_list(value)
                elif key != 'bins':
                    setattr(self, key, value)
            # set binaries
            if 'bins' in args:
                for key, value in args['bins'].items():
                    setattr(self.bins, key, value)
                self.bins = AttribDict(self.bins)

            # Workaround do reset bins
            if 'modes_cc' not in args and intype == 'CST':
                self.bins['mdcpl'] = 'matrixcst'
                self.bins['cst_partials'] = 'compu-deriv-cst'

        # Convert self.segmentsdir to AttribDict
        if len(self.segmentsdir) != 0:
            self = _get_data_dir(self)
            if sanity_checks is True:
                self.segmentsdir = check_segmentsdir(self)
                self.events, self.seg_suffix = get_events(self)
                self.model_size = _model_size(self)
        # if intype == 'CST':
        #     if not hasattr(self, 'mzero'):
        #         self.mzero = buildmodel(self)

    def __str__(self, extended=False):
        _pstr = "Setup for "
        if self.intype is not None:
            _pstr += "%s " % self.intype
        _pstr += "Inversion \n\n"

        _pstr += self.info()
        return _pstr

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def get_modes_in(self):
        modes_in = [str(len(self.modes))]
        for m, s in self.modes_sc.items():
            dnd = split_digit_nondigit(m)
            mode = ''
            for entry in dnd:
                entry = entry.lstrip('0')
                if len(entry) == 0:
                    entry = '0'
                mode += " %s" % entry
            if not mode.endswith('0') and s == 0:
                continue
            mode += ' %i 0' % s
            modes_in += [mode.upper()]

        if len(self.modes_cc) > 0 and len(self.modes_sc_dst) > 0:
            raise IOError("You cannot use dst and cstC inversion together")

        if len(self.modes_cc) > 0:
            modes_ccin = [str(len(self.modes_cc))]
            for m, s in self.modes_cc.items():
                mcc = split_digit_nondigit(' '.join(m.split('-')))
                # Remove spaces
                mcc = [x.strip(' ') for x in mcc]
                mcc = list(filter(None, mcc))
                for i, x in enumerate(mcc):
                    x = x.lstrip('0')
                    if len(x) == 0:
                        x = '0'
                    mcc[i] = x
                mode = ' '.join(mcc)
                mode += ' %i %i' % (s[0], s[1])
                modes_ccin += [mode.upper()]
        else:
            modes_ccin = None

        return modes_in, modes_ccin

    def get_modes_in_dst(self):
        if len(self.modes_sc_dst) > 0:
            modes_scin_dst = [str(len(self.modes_sc))]
            for m, s in self.modes_sc_dst.items():
                dnd = split_digit_nondigit(m)
                mode = ''
                for entry in dnd:
                    entry = entry.lstrip('0')
                    if len(entry) == 0:
                        entry = '0'
                    mode += " %s" % entry
                mode += ' %i 0' % s
                modes_scin_dst += [mode.upper()]
        else:
            modes_scin_dst = None

        return modes_scin_dst

    def build_pbs_script(self, path):
        pass

    def build_bash_script(self, path):
        pass

    def write(self, path=None, format='json'):
        if path is None:
            filename = os.path.join(self.rundir, 'setup.%s' % format)
        else:
            filename = os.path.join(path, 'setup.%s' % format)

        if format == 'pickle':
            write_pickle(self, filename)

        if format == 'json':
            data = {}
            for key, values in inspect.getmembers(self):
                if key.startswith('_'):
                    continue
                # Check if it is a bound method, if so ignore
                if hasattr(getattr(self, key), '__self__'):
                    continue
                if key == 'mzero':
                    data[key] = {}
                    for coeff, mode in values.items():
                        data[key][coeff] = {}
                        for m, c in mode.items():
                            if type(c) is not list:
                                data[key][coeff][m] = c.tolist()
                            else:
                                data[key][coeff][m] = c

                elif key == 'modes':
                    continue

                elif key in ('modes_sc', 'modes_cc', 'events', 'seg_suffix',
                             'segmentsdir'):
                    if key == 'modes_sc':
                        data['modes'] = OrderedDict(values)

                    data[key] = OrderedDict(values)

                else:
                    if type(values) not in (dict, AttribDict):
                        data[key] = values
                    else:
                        data[key] = dict(values)

            data = {'setup': data}
            with open(filename, 'w') as fh:
                json.dump(data, fh)

    def info(self):
        """
        prints the dictionary containing all setup attributes and values.

        Maybe sort this a bit?
        """
        content = vars(self)
        # keys = content.keys()
        # keys.sort()

        msg = ''
        for key in content.keys():
            value = content[key]
            if key == 'events':
                msg += 'events:\n'
                for set, cmts in value.items():
                    msg += '  %s: %i\n' % (set, len(cmts))
            elif key == 'mzero':
                msg += 'mzero:\n'
                if value is None:
                    msg += '   No Startupmodel defined\n'
                else:
                    for kind, modes in value.items():
                        msg += '  %s:\n' % kind
                        if modes is None:
                            msg += '   No Startupmodel defined\n'
                            continue
                        for mode, czero in modes.items():
                            mval = '  %s' % czero[:4]
                            msg += '    %s %s ...\n' % (mode, mval[:-1])
            elif key == 'segmentsdir':
                msg += 'segmentsdir:\n'
                for set, folder in value.items():
                    msg += '  %s: %s\n' % (set, folder)
            else:
                msg += '%s: %s\n' % (key, value)
        return msg

    def _internal_add_processing_info(self, info):
        """
        Add the given informational string to the `processing`
        """
        proc = self.processing
        proc.append(info)


def _get_data_dir(setup):
    if setup.rundir.startswith('/quanta'):
        setup.datadir = setup.default_data['quanta']

    elif setup.rundir.startswith('/scratch'):
        setup.datadir = setup.default_data['scratch']

        if hasattr(setup, 'bin_suffix') is True:
            for b in ['synseis_inv', 'cst_partials']:
                if not setup.bins[b].endswith('scratch'):
                    setup.bins[b] += '-scratch'

    elif setup.rundir.startswith('/nfs'):
        setup.datadir = setup.default_data['stig']
    else:
        setup.datadir = setup.default_data['local']
    pwd = os.getcwd()
    if (
        pwd.startswith('/quanta') or
        pwd.startswith('/scratch')
       ):
        if not os.path.exists(setup.datadir):
            msg = '{data} not found'.format(data=setup.datadir)
            raise IOError(msg)
    return setup


def _model_size(setup):
    if setup.intype == 'CST':
        N = 0
        for mode, smax in setup.modes_sc.items():
            ldeg = int(split_digit_nondigit(mode)[-1])
            # radial modes and other modes with smax > 0
            if ((smax > 0 and ldeg > 0) or (ldeg == 0)):
                N += N_splitting_coeffs(smax, 0)
                # Add d00
                N += 1

        for mode, s in setup.modes_cc.items():
            smin = s[0]
            smax = s[1]

            # CC with smax > 0
            if smax > 0:
                N += N_splitting_coeffs(smax, smin)

        for mode, smax in setup.modes_sc_dst.items():
            # dst always begins at 2
            if smax >= 2:
                N += N_splitting_coeffs(smax, 2)

    elif setup.intype == 'FSI':
        N = setup.model_size

    return N


def get_value_list(value):
    if type(value) in [str, int, float]:
        value = [value]
    return value


def get_mode_list(x):
    """
    x: dictionary like
    """
    mlist = []
    for mode, sdeg in x.items():
        mlist.append(mode)
    return mlist


def check_segmentsdir(setup):
    if (
      isinstance(setup.segmentsdir, AttribDict) or
      isinstance(setup.segmentsdir, dict)
         ):
        # add check for existing dirs here?
        for key, value in setup.segmentsdir.items():
            if not os.path.exists(value):
                msg = 'Segments not found: %s' % value
                raise IOError(msg)
        return setup.segmentsdir

    # If folder
    elif os.path.isdir(setup.segmentsdir):
        segname = basename(setup.segmentsdir)
        sdir = setup.segmentsdir
        return AttribDict({segname: sdir})

    # If file
    elif os.path.exists(setup.segmentsdir):
        segname = basename(setup.segmentsdir).split('.')[0]
        sdir = setup.segmentsdir
        return AttribDict({segname: sdir})

    else:
        msg = 'Segments not found: %s' % setup.segmentsdir
        raise IOError(msg)


def get_events(setup):

    accepted_ftypes = ['segmentZ', 'segmentT', 'segmentR', 'dat']

    if setup.rundir.startswith('/quanta'):
        accepted_cmts = os.listdir(setup.default_data['quanta'])
    elif setup.rundir.startswith('/scratch'):
        accepted_cmts = os.listdir(setup.default_data['scratch'])
    elif setup.rundir.startswith('/stig'):
        accepted_cmts = os.listdir(setup.default_data['stig'])
    else:
        accepted_cmts = os.listdir(setup.default_data['local'])

    suffix = AttribDict()
    segments = AttribDict()

    if type(setup.segmentsdir) in (dict, OrderedDict, AttribDict):
        for segname, sdir in setup.segmentsdir.items():
            segs = glob.glob('%s/[0-9]*[A-Z].*' % sdir)
            segs2 = []
            for seg in segs:
                if os.stat(seg).st_size != 0:
                    cmt, ftype = basename(seg).split('.')
                    if ftype not in accepted_ftypes:
                        continue
                    if "%s.cmt" % cmt not in accepted_cmts:
                        continue
                    if setup.use_Qfiles is False and cmt.endswith('Q'):
                        continue
                    if setup.use_Qfiles is True:
                        for _i, _s in enumerate(segs2):
                            if cmt[:6] == _s[:6]:
                                if cmt.endswith('Q'):
                                    segs2[_i] = cmt
                                if _s.endswith('Q'):
                                    continue

                    if cmt in segs2:
                        msg = 'Segment with multiple suffixes found:\n'
                        msg += '%s' % seg
                        raise IOError(msg)
                    segs2.append(cmt)

            segments[segname] = segs2
            suffix[segname] = ftype

    return segments, suffix


def modes_sanity_check(key, value_dict):
    if key == 'modes':
        for k, v in value_dict.items():
            if type(v) not in (int, np.int64):
                print(key, k, v)
                raise IOError('modes must have integer smax')
    elif key == 'modes_cc':
        for k, v in value_dict.items():
            if type(v) is not list:
                print(key, v)
                raise IOError('modes_cc must have [smin, smax]')
            _k = split_digit_nondigit(k)
            m = _k[:3] + _k[4:]

            if v[0] == 0 and v[-1] == 0:
                continue
            for _v in v:
                if _v not in max_cc_degrees(m):
                    print(k, v)
                    msg = 'mode pair couples only for these degrees: %s'
                    msg = msg % max_cc_degrees(m)
                    raise IOError(msg)
    elif key == 'modes_sc_dst':
        for k, v in value_dict.items():
            if type(v) is not int or v < 2:
                print(key, v)
                raise IOError('modes must have integer smax > 2')
        if len(value_dict) >= 2:
            raise IOError('Right now: dst inversion can only handle one mode')


def update_setup_pickle(setup):
    input = {'modes': {'0T2': 2}, 'startmodel': setup.model}
    setup_new = Setup(setup.intype, input)
    try:
        if 'cst' not in setup.mzero:
            setup_new.mzero['cst'] = setup.mzero
        else:
            setup_new.mzero = setup.mzero
    except AttributeError:
        setup_new.mzero = None

    for key, value in inspect.getmembers(setup):
        if key == 'mzero' or key.startswith('_'):
            continue
        if hasattr(getattr(setup, key), '__self__'):
            continue
        setattr(setup_new, key, value)

    return setup_new


def read(filename, format=None):
    # Try to figure out the file format
    if format is None:
        found = False
        if not filename.endswith('.json') and not filename.endswith('.pickle'):
            if os.path.exists(filename + '.json'):
                filename = filename + '.json'
                found = True
            if found is False:
                if os.path.exists(filename + '.pickle'):
                    filename = filename + '.pickle'
                else:
                    msg = 'No Setup file found'
                    raise IOError(msg)

        if filename.endswith('pickle'):
            format = 'pickle'

        elif filename.endswith('json'):
            format = 'json'
        else:
            msg = 'No Setup file found'
            raise IOError(msg)

    if format == 'pickle':
        return read_pickle(filename)

    elif format == 'json':
        with open(filename) as json_file:
            data = json.load(json_file, object_pairs_hook=OrderedDict)

        intype = data['setup']['intype']
        data['setup'].pop('intype')
        try:
            data['setup']['bins'] = AttribDict(data['setup']['bins'])
        except KeyError:
            pass
        try:
            uncertainties = AttribDict(data['setup']['uncertainties'])
            data['setup']['uncertainties'] = uncertainties
        except KeyError:
            pass
        setup = Setup(intype, args=data['setup'], sanity_checks=False)

    return setup


def update_setup(setup, host='scratch', update_bins=False):
    # if host == 'scratch':
    #     setup.rundir = setup.rundir.replace('quanta1/home', 'scratch')
    #     setup.datadir = setup.default_data['quanta']
    #     if update_bins is True:
    #         for b in ['synseis_inv', 'cst_partials']:
    #             if setup.bins[b].endswith('-quanta'):
    #                 setup.bins[b] = setup.bins[b].replace('-quanta',
    #                                                       '-scratch')
    #             if not setup.bins[b].endswith('-scratch'):
    #                 setup.bins[b] = setup.bins[b] + '-scratch'
    #
    # elif host == 'quanta':
    #     setup.rundir = setup.rundir.replace('scratch', 'quanta1/home')
    #     if update_bins is True:
    #         for b in ['synseis_inv', 'cst_partials']:
    #             if setup.bins[b].endswith('quanta'):
    #                 setup.bins[b] = setup.bins[b].replace('-quanta', '')

    if len(setup.modes_cc) == 0 and len(setup.modes_sc) > 1:
        modes = []
        for m, deg in setup.modes_sc.items():
            if deg != 0:
                modes = [(m, deg)] + modes
            else:
                modes.append((m, deg))
        setup.modes_sc = OrderedDict(modes)

    return setup
