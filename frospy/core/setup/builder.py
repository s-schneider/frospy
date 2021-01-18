#!/usr/bin/env python2
from obspy.core import AttribDict
from frospy.core.setup.settings import Setup, update_setup
from frospy.core.setup.settings import read as read_setup
from frospy.core.setup import create
import shutil
import glob
from os.path import exists, join

python_path = "/home/simons/dev/nmPy/nmpy/standalone_tools"
merge_python = "%s/merge.py" % python_path


def build_synthetics(args=None, verbose=True, remove_existing=False,
                     run_code=False):
    """
    setup for calculating synthetics.
    """
    setup = Setup('synthetics')
    if remove_existing is True:
        if exists(setup.rundir):
            shutil.rmtree(setup.rundir)

    # Setup config folder ?
    # Setup data, cmt, segments, modes.in etc.
    # Setup bash and pbs scripts

    setup = AttribDict()
    return setup


def build_inversion(input_dict=None, inversion_type='CST',
                    remove_existing=False, submit=False, output='pbs',
                    check_segments=False):
    """
    input_dict:

    Needed input_dict keys:
     - rundir
     - startmodel
     - modes
     - damping
     - iterations

    Optional input_dict keys:

     - bindir
     - bins
     - channeldamping
     - cmtdir
     - cross_coupling
     - damping
     - datadir
     - ellip
     - intype
     - iterno
     - model
     - pbsjob
     - rotation
     - run_keys
     - segmentsdir
     - s_deg
     - splittingdir     # FSI
     - modedirs         # FSI

    Example:
    from frospy.core.setup.builder import build_inversion
    from obspy.core import AttribDict
    from collections import OrderedDict

    input = AttribDict()
    mdir = '//nfs/stig/simons/splitting/modes/00t04/'
    input['rundir'] = os.path.join(mdir, 'self_coupling/Runs/T-Comp/invtest')
    input['bindir'] = '/home/simons/bin'
    input['startmodel'] = 'PREM'
    input['comp'] = 'T'
    input['segmentsdir'] = os.path.join(mdir, 'segments_T')
    input['damping'] = np.logspace(-3, 1, 5).tolist()
    input['iterations'] = [1, 10]
    input['modes'] = OrderedDict([("0T4", 4)])
    input['pbsjob'] ='invtest'


    setup = build_inversion(input, 'CST', remove_existing=True,
                            seg_suffix=None)

    """

    # Create setup objects with all defaults
    sanity_check_input(input_dict)
    setup = Setup(inversion_type, input_dict)
    sanity_check_setup(setup, check_segments)

    # create rundir
    create.rundir(setup, remove_existing)

    # create modes.in
    if setup.intype == 'CST':
        create.modesin(setup)

    # create segment_files
    create.segments(setup)

    # create startmodel
    create.startmodel(setup)

    # create script
    create.script(setup, output)

    # submit job
    if submit is True:
        create.submission(setup)

    setup.write()
    return setup


def build_array_inversion(input_dict=None, inversion_type='CST',
                          remove_existing=False, check_segments=False,
                          scripts_only=False, iter_damp_test=False):
    """
    input_dict:

    Needed input_dict keys:
     - rundir
     - startmodel
     - modes
     - damping
     - iterations

    Optional input_dict keys:

     - bindir
     - bins
     - channeldamping
     - cmtdir
     - cross_coupling
     - damping
     - datadir
     - ellip
     - intype
     - iterno
     - model
     - pbsjob
     - rotation
     - run_keys
     - segmentsdir
     - s_deg
     - splittingdir     # FSI
     - modedirs         # FSI

    Example:
    from frospy.core.setup.builder import build_array_inversion
    from obspy.core import AttribDict

    input = AttribDict()
    mdir = '//nfs/stig/simons/splitting/modes/00t04/'
    input['rundir'] = os.path.join(mdir, 'self_coupling/Runs/T-Comp/test')
    input['bindir'] = '/home/simons/bin'
    input['startmodel'] = 'PREM'
    input['comp'] = 'T'
    input['segmentsdir'] = os.path.join(mdir, 'segments_T')
    input['damping'] = np.logspace(-3, 1, 5).tolist()
    input['iterations'] = [1, 10]
    input['modes'] = {'0T4': 4}

    setup = build_array_inversion(input, 'CST', remove_existing=True)

    """

    if type(input_dict) != Setup:
        sanity_check_input(input_dict)
        setup = Setup(inversion_type, input_dict)
    else:
        setup = input_dict

    sanity_check_setup(setup, check_segments)

    if scripts_only is not True:
        # create rundir
        create.rundir(setup, remove_existing)

        # create modes.in
        if setup.intype == 'CST':
            create.modesin(setup)

        # create eventsdir and segment_files
        create.arraydirs(setup, remove_existing)

        # create startmodel
        create.startmodel(setup)

    # create the 3 pbs_scripts
    create.startup_script(setup)
    create.array_script(setup)
    create.new_script(setup)

    # create dependencies script
    create.dependencies_script(setup, iter_damp_test=iter_damp_test)
    if setup.intype == 'CST':
        setup.write()
    return setup


def build_uncertainties(rundir, damping, remove_existing=False,
                        size_of_subsets='full', N_of_subsets=200,
                        allevents_subset_ratio=None,
                        jack_knive=True,
                        boot_strap=False):
    """
    Test script in:
    /scratch/simons/splitting/modes/00t04/self_coupling/TR-Comp/deg8

    sbatch submit_uncertainties.sh

    """
    setup = read_setup(join(rundir, 'setup'))

    setup = create.uncertainties_dir(setup, remove_existing,
                                     size_of_subsets=size_of_subsets,
                                     N_of_subsets=N_of_subsets,
                                     allevents_subset_ratio=allevents_subset_ratio,
                                     jack_knive=jack_knive,
                                     boot_strap=boot_strap)
    setup = create.uncertainties_script(setup, damping)
    setup.write()

    return setup


def sanity_check_input(input_dict):
    checkpoints = ['rundir', 'startmodel', 'damping',
                   'iterations', 'modes']

    for point in checkpoints:
        if point not in input_dict:
            raise IOError("%s not defined in input_dict" % point)
    return


def sanity_check_setup(setup, check_segments=False):
    if setup.modes is None:
        raise IOError('No modes specified!')

    if setup.startmodel is None:
        raise IOError('No inmodel specified')

    if check_segments is True:
        N = 0
        for segdir in setup.segmentsdir.values():
            N += len(glob.glob("%s/*.dat" % segdir))
            N += len(glob.glob("%s/*.segment*" % segdir))
        if N == 0:
            raise IOError('No segments Found')
    return


def update_scripts(rundir, host='scratch', intype='CST'):
    setup = read_setup(join(rundir, 'setup.json'))
    setup = update_setup(setup, host)

    for b in ('synseis_inv', 'cst_partials'):
        if setup.bins[b].endswith('-scratch') and host == 'scratch':
            setup.bins[b] = setup.bins[b].replace('-scratch', '')
    if setup.bins[b].endswith('-quanta') and host == 'scratch':
        setup.bins[b] = setup.bins[b].replace('-quanta', '')

    create.rundir(setup, remove_existing=False, update_only=True)
    setup = build_array_inversion(setup, intype, remove_existing=False,
                                  scripts_only=True)
    return setup
