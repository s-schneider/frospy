import os
from os.path import exists, join, splitext
import shutil
from shutil import copy as cp
import subprocess
import numpy as np
from obspy.core import AttribDict

from frospy.core.modes import read as read_modes
from frospy.core.modes import format_name
from frospy.core.segment import read as read_seg

import frospy.standalone_tools as pypath


invoutpath = join(pypath.__path__[0], 'write_inversion_out.py')
if invoutpath.startswith('/net'):
    invoutpath = 'nmpy/standalone_tools/write_inversion_out.py'

mergepath = join(pypath.__path__[0], 'merge.py')
if mergepath.startswith('/net'):
    mergepath = 'nmpy/standalone_tools/merge.py'

forkpath = join(pypath.__path__[0], 'read_event_loop.pl')
if forkpath.startswith('/net'):
    forkpath = 'nmpy/standalone_tools/read_event_loop.pl'


def rundir(setup, remove_existing, update_only=False):
    """
    Creates the run directory as specifed in setup.rundir
    """
    if not update_only:
        if exists(setup.rundir):
            if remove_existing is True:
                shutil.rmtree(setup.rundir)
            else:
                raise IOError('rundir exists!')
        os.makedirs(setup.rundir)

        # If inversion_out is default, then it is within rundir
        if setup.inversion_outdir == 'inversion_out':
            os.makedirs(join(setup.rundir, setup.inversion_outdir))
        else:
            os.makedirs(setup.inversion_outdir)

    # Copy write_inversion_out.py to folder
    src = [join(setup.nmpydir, invoutpath),
           join(setup.nmpydir, mergepath),
           join(setup.nmpydir, forkpath)]
    dst = [join(setup.rundir, "write_inversion_out.py"),
           join(setup.rundir, "merge.py"),
           join(setup.rundir, "read_event_loop.pl")]

    for s, d in zip(src, dst):
        shutil.copy(s, d)

    return


def arraydirs(setup, remove_existing):
    """
    Creates a directory for every event specifed in setup.segmentsdir
    with the modes.in and segment file
    """
    if (
        isinstance(setup.segmentsdir, AttribDict) or
        isinstance(setup.segmentsdir, dict)
       ):
        nof_sets = len(setup.segmentsdir)
        if setup.intype == 'CST':
            fname = os.path.join(setup.rundir, 'allevents')

            with open(fname, 'w+') as fh:
                for segname, sdir in setup.segmentsdir.items():
                    events = sorted(setup.events[segname])

                    for event in events:
                        if nof_sets > 1:

                            segdir = os.path.join(setup.rundir, segname)
                            if not exists(segdir):
                                os.makedirs(segdir)
                            eventdir = os.path.join(setup.rundir, segname,
                                                    event)

                        else:
                            eventdir = os.path.join(setup.rundir, event)

                        if exists(eventdir):
                            if remove_existing is True:
                                shutil.rmtree(eventdir)
                            else:
                                raise IOError('eventdir exists!')
                        os.makedirs(eventdir)
                        # write allevents
                        if nof_sets > 1:
                            fh.write('%s %s\n' % (segname, event))
                        else:
                            fh.write('%s\n' % event)

                        # copying modesin
                        modesfile = modesin(setup, filename='modes.in')
                        cp(modesfile, eventdir)

            # copying segments from rundir to each eventdir
            src_all, dst_all = segments(setup, print_seg=False)
            for src, dst in zip(src_all, dst_all):
                dst = splitext(dst)[0]
                cp(src, dst)

        elif setup.intype == 'FSI':
            fname = os.path.join(setup.rundir, 'allevents')
            with open(fname, 'w+') as fh:
                for segname, sdir in setup.segmentsdir.items():
                    events = sorted(setup.events[segname])
                    for event in events:
                        eventdir = os.path.join(sdir, event)
                        if exists(eventdir):
                            if remove_existing is True:
                                shutil.rmtree(eventdir)
                            else:
                                raise IOError('eventdir exists!')
                        os.makedirs(eventdir)
                        # write allevents
                        fh.write('%s %s\n' % (segname, event))

                        # copying modesin
                        modesfile = os.path.join(sdir, 'modes.in')
                        modessfile = os.path.join(sdir, 'modess.in')
                        modesmfile = os.path.join(sdir, 'modesm.in')
                        segfile = os.path.join(sdir, '%s.dat' % (event))
                        cp(modesfile, eventdir)
                        cp(modessfile, eventdir)
                        cp(modesmfile, eventdir+'/modes.in')
                        cp(segfile, eventdir)

    elif os.path.isdir(setup.segmentsdir):
        fname = os.path.join(setup.rundir, 'allevents')
        events = sorted(setup.events)
        with open(fname, 'w+') as fh:
            for event in events:
                eventdir = os.path.join(setup.rundir, event)
                if exists(eventdir):
                    if remove_existing is True:
                        shutil.rmtree(eventdir)
                    else:
                        raise IOError('eventdir exists!')
                os.makedirs(eventdir)
                fh.write('%s\n' % event)  # write allevents

                # copying modesin
                modesfile = modesin(setup, filename='modes.in')
                cp(modesfile, eventdir)

        # copying segments from rundir to each eventdir
        src_all, dst_all = segments(setup, print_seg=False)
        for src, dst in zip(src_all, dst_all):
            dst = splitext(dst)[0]
            cp(src, dst)

    return


def modesin(setup, filename='modes.in'):
    """
    Writes modes.in file needed for mdcplmrho* programs
    """
    path = setup.rundir
    fname = os.path.join(path, filename)
    with open(fname, 'w+') as fh:
        fh.write('%s\n' % len(setup.modes_sc))
        for mode, smax in setup.modes_sc.items():
            m = read_modes(modenames=mode)[0]
            fh.write("%03d %s %03d\n" % (m.n, m.type.lower(), m.l))
    return fname


def segments(setup, print_seg=True):
    """
    Setup of segment files. Copies all files ending with 'suffix'
    from setup.segmentsdir to setup.rundir .
    """
    dst_all = []
    src_all = []
    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch')
       ):
        default = '/quanta1/home/simons/splitting/alldatafiles'
    elif setup.rundir.startswith('/nfs'):
        default = '//nfs/stig/deuss/modes/alldatafiles'
    else:
        default = setup.datadir


    accepted_ftypes = ['segmentZ', 'segmentT', 'segmentR', 'dat']
    if setup.segmentsdir is None:
        return

    # If segmentsdir is a dict-type
    if (
        isinstance(setup.segmentsdir, AttribDict) or
        isinstance(setup.segmentsdir, dict)
       ):
        nof_sets = len(setup.segmentsdir)
        for segname, sdir in setup.segmentsdir.items():
            for seg in os.listdir(sdir):
                cmt, ftype = seg.split('.')
                if ftype not in accepted_ftypes:
                    continue

                if "%s.cmt" % cmt not in os.listdir(default):
                    continue
                src = join(sdir, seg)
                if nof_sets > 1:
                    dst = join(setup.rundir, segname, seg)
                else:
                    dst = join(setup.rundir, seg)
                dst_all.append(dst)
                src_all.append(src)
                if not exists(dst) and print_seg is True:
                    try:
                        cp(src, dst)
                    except Exception:
                        print("Needs array inversion setup")

    # If segmentsdir is a directory
    elif os.path.isdir(setup.segmentsdir):
        for seg in os.listdir(setup.segmentsdir):
            cmt, ftype = seg.split('.')
            if ftype not in accepted_ftypes:
                continue

            if "%s.cmt" % cmt not in os.listdir(default):
                continue
            src = join(setup.segmentsdir, seg)
            dst = join(setup.rundir, seg)
            dst_all.append(dst)
            src_all.append(src)
            if not exists(dst) and print_seg is True:
                cp(src, dst)

    # If segmentsdir is a file
    else:
        src = join(setup.segmentsdir)
        seg = setup.segmentsdir.split('/')[-1]
        cmt, ftype = seg.split('.')
        if ftype in accepted_ftypes and print_seg is True:
            if "%s.cmt" % cmt in os.listdir(default):
                dst = join(setup.rundir, seg)
                dst_all.append(dst)
                src_all.append(src)
                cp(src, dst)

    return src_all, dst_all


def submission(setup):
    """
    Submits the pbs script to stig
    """
    pbsrundir = pbsrundir_sanity_check(setup.rundir)
    msg = "ssh stig 'cd %s; qsub qsub_me.pbs'" % (pbsrundir)
    p = subprocess.Popen(msg, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    if err is not None:
        raise RuntimeError("Submission Failed:\n{}".format(output))


def damping_file(setup, mdamp):
    modelfile = 'Cm1.dat'

    if setup.intype == 'CST':
        mfile = join(setup.rundir, modelfile)
        with open(mfile, 'w+') as fh:
            for key, val in mdamp['cst'].items():
                for v in val:
                    fh.write('%s\n' % v)

            if len(setup.modes_sc_dst) > 0:
                for key, val in mdamp['dst'].items():
                    for v in val:
                        fh.write('%s\n' % v)
    return


def startmodel(setup):
    """
    Setup of startmodel.
    creates the startfiles file in setup.rundir, reading from setup.mzero
    """
    if setup.intype == 'CST':
        modelfile = 'mzero.dat'

    if setup.startmodel == 'S20RTS':
        if setup.intype == 'CST':  # PREM & S20RTS are the same for cst
            mfile = join(setup.rundir, modelfile)
            with open(mfile, 'w+') as fh:
                for key, val in setup.mzero["cst"].items():
                    for v in val:
                        fh.write('%s\n' % v)

                if len(setup.modes_sc_dst) > 0:
                    for key, val in setup.mzero["dst"].items():
                        for v in val:
                            fh.write('%s\n' % v)
            # To start, mcst.dat must be there too
            cp(mfile, join(setup.rundir, 'mcst.dat'))
            setup.model = mfile

        elif setup.intype == 'FSI':
            if setup.s_deg == 20:
                if (
                    setup.rundir.startswith('/quanta') or
                    setup.rundir.startswith('/scratch')
                   ):
                    setup.model = '/quanta1/home/jagt/dta/S20RTS_MC.sph'
                else:
                    setup.model = '/home/deuss/dta/S20RTS_MC.sph'
            elif setup.s_deg == 4:
                if (
                    setup.rundir.startswith('/quanta') or
                    setup.rundir.startswith('/scratch')
                   ):
                    setup.model = '/quanta1/home/jagt/dta/S20RTS_l4.sph'
                else:
                    setup.model = '/home/deuss/dta/S20RTS_l4.sph'
            else:
                pass

    elif setup.startmodel == 'PREM':
        if setup.intype == 'CST':
            mfile = join(setup.rundir, modelfile)
            with open(mfile, 'w+') as fh:
                for key, val in setup.mzero["cst"].items():
                    for v in val:
                        fh.write('%s\n' % v)
                if len(setup.modes_sc_dst) > 0:
                    for key, val in setup.mzero["dst"].items():
                        for v in val:
                            fh.write('%s\n' % v)
            # To start from PREM, mcst.dat must be there too
            cp(mfile, join(setup.rundir, 'mcst.dat'))
            setup.model = mfile

        elif setup.intype == 'FSI':
            if setup.s_deg == 20:
                if (
                    setup.rundir.startswith('/quanta') or
                    setup.rundir.startswith('/scratch')
                   ):
                    setup.model = '/quanta1/home/jagt/dta/mzero_l20.sph'
                else:
                    setup.model = '/home/jagt/splitting/mzero_l20.sph'
            elif setup.s_deg == 4:
                if (
                    setup.rundir.startswith('/quanta') or
                    setup.rundir.startswith('/scratch')
                   ):
                    setup.model = '/quanta1/home/jagt/dta/mzero_l4.sph'
                else:
                    setup.model = '/home/jagt/splitting/mzero_l4.sph'
            else:
                pass
    else:
        if setup.intype == 'CST':
            mfile = join(setup.rundir, modelfile)
            with open(mfile, 'w+') as fh:
                for key, val in setup.mzero["cst"].items():
                    for v in val:
                        fh.write('%s\n' % v)
                if len(setup.modes_sc_dst) > 0:
                    for key, val in setup.mzero["dst"].items():
                        for v in val:
                            fh.write('%s\n' % v)
            # To start from PREM, mcst.dat must be there too
            cp(mfile, join(setup.rundir, 'mcst.dat'))
            setup.model = mfile

    return


def get_damping_string(setup):
    damping = ''
    for d in setup.damping:
        damping += "%4.2e " % d
    return damping


def get_segments(seg_suffix):
    msg = ''
    msg += 'declare -a SEGMENTS='
    msg += '(`ls | egrep ^[0-9]*[A-Z].%s`)\n' % seg_suffix
    return msg


def get_FILE_extension(seg_suffix):
    msg = "${{FILE%%%%.%s}}" % seg_suffix
    return msg


def dorun_preprundir(setup):
    msg = 'cp {zeromodel} .\n'
    if setup.model.endswith('sph'):
        msg += '{sph2m} <<EOF > mzero.dat\n'
        msg += '{zeromdl}\n'
        msg += 'EOF\n'
    else:
        msg += 'cp {zeromodel} mzero.dat\n'
    msg += 'cp mzero.dat mold.dat\n'
    msg += '\n'
    msg = msg.format(sph2m=join(setup.bindir, setup.bins.sph2m),
                     zeromodel=setup.model,
                     zeromdl=setup.model.split('/')[-1])
    return msg


def dorun_startup(setup):
    # Startmodel mzero.dat and mcst.dat are generated by create.startmodel
    msg = '#---------------------------#\n'
    msg += '#  STARTUP WITH ZERO MODEL  #\n'
    msg += '#---------------------------#\n'
    msg += 'echo "Starting up" >> progress.out\n'
    msg += 'rm inversion.dat inversion.out 2> /dev/null\n'
    msg += 'rm dsyni.dat dsynr.dat dobsi.dat dobsr.dat mode.dat att.dat'
    msg += 'misfit.dat 2> /dev/null\n'

    if setup.intype == 'CST':
        m = 'cp mzero.dat mcst.dat\n\n'
        m += '{rot_ellip}\n'
        m = m.format(rot_ellip=join(setup.bindir, setup.bins.rot_ellip))

        if len(setup.modes_sc_dst) > 0:
            m += '{write_dst} << EOF\n'
            for mode, smax in setup.modes_sc.items():
                m += '%s\n' % smax
            for mode, smax in setup.modes_sc_dst.items():
                m += '%s\n' % smax
            m += 'EOF\n\n'
            m = m.format(write_dst=join(setup.bindir, setup.bins.write_dst))

        m += '{mdcpl} << EOF\n'
        for mode, smax in setup.modes_sc.items():
            # /home/deuss/bin/matrixcst
            # maximum structural deg for every mode
            m += '0 %s\n' % smax
        if hasattr(setup, 'modes_cc'):
            # /home/deuss/bin/matrixcstC
            for mode, s in setup.modes_cc.items():
                m += '%s %s\n' % (s[0], s[1])
        m += 'EOF\n'
        m = m.format(mdcpl=join(setup.bindir, setup.bins.mdcpl))

        if len(setup.modes_sc_dst) > 0:
            m += 'mv matrix.dat matrix-cst.dat\n\n'
            # /home/deuss/bin/matrixdst
            m += '{mdcplmq} << EOF\n'
            for mode, smax in setup.modes_sc_dst.items():
                # /home/deuss/bin/matrixcst
                # maximum structural deg for every mode
                m += '2 %s\n' % smax
            m += 'EOF\n'
            m += 'mv matrix.dat matrix-dst.dat\n\n'
            m += '{addmdcplmq} << EOF\n'
            m += 'matrix-cst.dat\n'
            m += 'matrix-dst.dat\n'
            m += 'EOF\n'
            m = m.format(mdcplmq=join(setup.bindir, setup.bins.mdcplmq),
                         addmdcplmq=join(setup.bindir, setup.bins.addmdcplmq))

    elif setup.intype == 'FSI':
        m = '{mdcpl} -model {zeromodel} <<EOF\n'
        m += '0 {s_deg}\n'
        m += '{rotation}\n'
        m += '{ellip}\n'
        m += 'EOF\n'
        m = m.format(mdcpl=join(setup.bindir, setup.bins.mdcpl),
                     zeromodel=setup.model, s_deg=setup.s_deg,
                     rotation=setup.rotation, ellip=setup.ellip)

    msg += m
    msg += '\n'
    msg += '{matdiag} << EOF\n'  # /home/deuss/bin/matdiag
    msg += '0\n'
    msg += 'EOF\n\n'

    msg += get_segments(list(setup.seg_suffix.values())[0])
    msg += 'for FILE in "${{SEGMENTS[@]}}"\n'
    msg += 'do\n'
    msg += '{synseis_inv} '
    msg += get_FILE_extension(list(setup.seg_suffix.values())[0])
    msg += ' << EOF\n'
    msg += 'rr.dat\n'
    msg += 'rr1.dat\n'
    msg += 'ww.dat\n'
    msg += 'EOF\n'
    msg += 'done\n\n'
    if setup.intype == 'CST':
        msg += 'cp mcst.dat mold.dat\n'
        msg += 'cp mzero.dat {summary_dir}/mzero.dat\n'
    msg += '{avmisfit}\n'  # home/deuss/bin/avmisfit
    msg += 'python {invout} $RUNDIR 0 {summary_dir} {schema}\n'

    msg = msg.format(matdiag=join(setup.bindir, setup.bins.matdiag),
                     synseis_inv=join(setup.bindir, setup.bins.synseis_inv),
                     avmisfit=join(setup.bindir, setup.bins.avmisfit),
                     invout="write_inversion_out.py",
                     summary_dir=setup.inversion_outdir,
                     schema=setup.schema)
    return msg


def dorun_cst_partials(setup):
    msg = '\n'
    msg += '#----------------#\n'
    msg += '#  Cst Partials  #\n'
    msg += '#----------------#\n'
    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'sed -i "/$iterno $damp $SEGDIR $EVENT/s/$/ '
            msg += 'partials/" $PO\n'
        else:
            msg += 'sed -i "/$iterno $damp $EVENT/s/$/ partials/" $PO\n'
    msg += 'rm Ar.dat Ai.dat allevents_partials 2> /dev/null\n'

    for mode, smax in setup.modes_sc.items():
        # /home/deuss/bin/compu-deriv-cst
        # maximum structural deg for every mode
        msg += 'echo %s >> allevents_partials\n' % smax

    if hasattr(setup, 'modes_cc'):
        # /home/deuss/bin/compu-deriv-cstC
        for mode, s in setup.modes_cc.items():
            msg += 'echo %s %s >> allevents_partials\n' % (s[0], s[1])

    msg += 'for FILE in "${{SEGMENTS[@]}}"\n'
    msg += 'do\n'
    msg += 'echo '
    msg += get_FILE_extension(list(setup.seg_suffix.values())[0])
    msg += ' >> allevents_partials\n'
    msg += 'done\n'
    msg += '{partial_derivatives} < allevents_partials\n'
    msg += 'count=`grep "NaN" Ai.dat | wc -l`\n'
    msg += '[ ! $count == 0 ] && echo '
    msg += '"    -NaN found in "'
    msg += '$SEGDIR/$EVENT " Ai.dat" >> progress.out && rm A* && '
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg = msg.format(partial_derivatives=join(setup.bindir,
                                              setup.bins.cst_partials))
    return msg


def dorun_dst_partials(setup):
    msg = '#----------------#\n'
    msg += '#  Dst Partials  #\n'
    msg += '#----------------#\n'
    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'sed -i "/$damp $SEGDIR $EVENT/s/$/ dst-partials/" $PO\n'
        else:
            msg += 'sed -i "/$damp $EVENT/s/$/ dst-partials/" $PO\n'
    msg += 'rm allevents_partials\n'

    for mode, smax in setup.modes_sc_dst.items():
        msg += 'echo %s >> allevents_partials\n' % smax

    msg += 'for FILE in "${{SEGMENTS[@]}}"\n'
    msg += 'do\n'
    msg += 'echo '
    msg += get_FILE_extension(list(setup.seg_suffix.values())[0])
    msg += ' >> allevents_partials\n'
    msg += 'done\n'
    msg += '{partial_derivatives} < allevents_partials\n'
    msg += 'count=`grep "NaN" Ai.dat | wc -l`\n'
    msg += '[ ! $count == 0 ] && echo '
    msg += '"    -NaN found in "'
    msg += '$SEGDIR/$EVENT " Ai.dat" >> progress.out && rm A* && '
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg = msg.format(partial_derivatives=join(setup.bindir,
                                              setup.bins.dst_partials))
    return msg


def dorun_cst_partials_array(setup):
    msg = '#----------------#\n'
    msg += '#  Cst Partials  #\n'
    msg += '#----------------#\n'

    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'sed -i "/$iterno $damp $SEGDIR $EVENT/s/$/ '
            msg += 'partials/" $PO\n'
        else:
            msg += 'sed -i "/$iterno $damp $EVENT/s/$/ partials/" $PO\n'

    msg += 'rm Ar.dat Ai.dat allevents_partials 2> /dev/null\n'

    for mode, smax in setup.modes_sc.items():
        # /home/deuss/bin/compu-deriv-cst
        # maximum structural deg for every mode
        msg += 'echo %s >> allevents_partials\n' % smax
    if hasattr(setup, 'modes_cc'):
        # /home/deuss/bin/compu-deriv-cstC
        for mode, s in setup.modes_cc.items():
            msg += 'echo %s %s >> allevents_partials\n' % (s[0], s[1])

    msg += 'echo $EVENT >> allevents_partials\n\n'
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg += 'count=`grep "NaN" Ai.dat | wc -l`\n'
    msg += '[ ! $count == 0 ] && echo '
    msg += '"    -NaN found in "'
    msg += '$SEGDIR/$EVENT " Ai.dat" >> $MAINDIR/progress.out && rm A* && '
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg = msg.format(partial_derivatives=join(setup.bindir,
                                              setup.bins.cst_partials))
    return msg


def dorun_dst_partials_array(setup):
    msg = '#----------------#\n'
    msg += '#  Dst Partials  #\n'
    msg += '#----------------#\n'

    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'sed -i "/$damp $SEGDIR $EVENT/s/$/ dst-partials/" $PO\n'
        else:
            msg += 'sed -i "/$damp $EVENT/s/$/ dst-partials/" $PO\n'

    msg += 'rm allevents_partials\n'

    for mode, smax in setup.modes_sc_dst.items():
        msg += 'echo %s >> allevents_partials\n' % smax

    msg += 'echo $EVENT >> allevents_partials\n\n'
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg += 'count=`grep "NaN" Ai.dat | wc -l`\n'
    msg += '[ ! $count == 0 ] && echo '
    msg += '"    -NaN found in "'
    msg += '$SEGDIR/$EVENT " Ai.dat" >> $MAINDIR/progress.out && rm A* && '
    msg += '{partial_derivatives} < allevents_partials\n\n'
    msg = msg.format(partial_derivatives=join(setup.bindir,
                                              setup.bins.dst_partials))
    return msg


def dorun_buildATA_array(setup):
    msg = '#--------------------#\n'
    msg += '#  Build ATA matrix  #\n'
    msg += '#--------------------#\n'

    msg += 'if [ ! -f $EVENT.dat ]; then\n'
    msg += "        suffix=`ls $EVENT.segment? | cut -d'.' -f2`\n"
    msg += 'else\n'
    msg += '        suffix=dat\n'
    msg += 'fi\n'

    msg += 'N_seg=`cat $EVENT.$suffix | wc -l`\n'
    msg += 'if [ ! $N_seg == 0 ]; then\n'
    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'sed -i "/$iterno $damp $SEGDIR $EVENT/s/$/ '
            msg += 'buildATA/" $PO\n'
        else:
            msg += 'sed -i "/$iterno $damp $EVENT/s/$/ buildATA/" $PO\n'

    msg += '{buildATA} << EOF\n'
    msg += '100.0\n'
    msg += 'EOF\n'
    msg += 'fi\n'
    msg = msg.format(buildATA=join(setup.bindir, setup.bins.buildATA))
    return msg


def dorun_fsi_partials(setup):
    msg = '\n'
    msg += '#----------------#\n'
    msg += '#  FSI Partials  #\n'
    msg += '#----------------#\n'
    msg += 'rm Ar.dat Ai.dat allevents\n'
    msg += get_segments(list(setup.seg_suffix.values())[0])
    msg += 'for FILE in "${{SEGMENTS[@]}}"\n'
    msg += 'do\n'
    msg += 'echo '
    msg += get_FILE_extension(list(setup.seg_suffix.values())[0])
    msg += ' >> allevents\n'
    msg += 'done\n'
    msg += '{partial_derivatives} < allevents\n'

    msg = msg.format(partial_derivatives=join(setup.bindir,
                                              setup.bins.fsi_partials))
    return msg


def dorun_fsi_partials_array(setup):
    msg = '#----------------#\n'
    msg += '#  Fsi Partials  #\n'
    msg += '#----------------#\n'

    msg += 'rm Ar.dat Ai.dat allevents\n'

    msg += 'echo $EVENT >> allevents\n\n'
    msg += '{partial_derivatives} < allevents\n'

    msg += '\n'
    msg += '#--------------------#\n'
    msg += '#  Build ATA matrix  #\n'
    msg += '#--------------------#\n'
    msg += '{buildATA} << EOF\n'
    msg += '1000.0\n'

    msg += 'EOF\n'

    msg = msg.format(buildATA=join(setup.bindir, setup.bins.buildATA),
                     partial_derivatives=join(setup.bindir,
                                              setup.bins.fsi_partials))
    return msg


def dorun_buildATA(setup):
    msg = '\n'
    msg += '#--------------------#\n'
    msg += '#  Build ATA matrix  #\n'
    msg += '#--------------------#\n'
    if setup.intype == 'CST':
        msg += 'echo "    -buildATA" >> progress.out\n'
    msg += '{buildATA} << EOF\n'
    if setup.intype == 'CST':
        msg += '100.0\n'
    elif setup.intype == 'FSI':
        msg += '1000.0\n'
    msg += 'EOF\n\n'
    msg = msg.format(buildATA=join(setup.bindir, setup.bins.buildATA))
    return msg


def dorun_addATA(setup, uncertainties=False):
    rundir = setup.rundir
    if rundir.startswith('//nfs'):
        rundir = rundir.replace('//nfs/stig/', '/home/')

    msg = '#--------------#\n'
    msg += '#  Adding ATA  #\n'
    msg += '#--------------#\n'
    if setup.intype == 'CST':  # for PBS array only
        if uncertainties is False:
            Ndir = '$RUNDIR/allevents'
        else:
            Ndir = '$set_name'

        # if len(setup.events) > 1:
        #     msg += 'echo "    -addATA" >> $RUNDIR/progress.out\n'
        #     msg += 'rm ATAmatrix.dat ATd.dat dTd.dat ATA.dat ATA.list\n\n'
        #     msg += 'Nmodes={msize}\n'
        #     msg += 'echo "$Nmodes" > ATA.list\n'
        #     msg += 'Nevents=`cat %s | wc -l`\n' % Ndir
        #     msg += 'echo "$(( Nevents-1 ))" >> ATA.list\n\n'
        #     msg += 'i=0\n'
        #     msg += 'for f in `cat %s | sed "s/ /\//g"`; do\n' % Ndir
        #     msg += 'if [ $i -eq 0 ]; then\n'
        #     msg += '    cp $RUNDIR/$f/ATAmatrix.dat ATA.dat\n'
        #     msg += '    cp $RUNDIR/$f/ATd.dat .\n'
        #     msg += '    cp $RUNDIR/$f/dTd.dat .\n'
        #     msg += 'else\n'
        #     msg += '    echo $RUNDIR/$f/ATAmatrix.dat >> ATA.list\n'
        #     msg += '    echo $RUNDIR/$f/ATd.dat >> ATA.list\n'
        #     msg += '    echo $RUNDIR/$f/dTd.dat >> ATA.list\n'
        #     msg += 'fi\n'
        #     msg += 'i=$(( i+1 ))\n'
        #     msg += 'done\n'
        #     msg += '{addATA} < ATA.list\n'
        #
        # else:
        #     msg += 'echo "    -addATA" >> $RUNDIR/progress.out\n'
        #     msg += 'rm ATAmatrix.dat ATd.dat dTd.dat ATA.dat ATA.list\n\n'
        #     msg += 'Nmodes={msize}\n'
        #     msg += 'echo "$Nmodes" > ATA.list\n'
        #     msg += 'Nevents=`cat %s | wc -l`\n' % Ndir
        #     msg += 'echo "$(( Nevents-1 ))" >> ATA.list\n\n'
        #     msg += 'i=0\n'
        #     msg += 'for f in `cat %s | sed "s/ /\//g"`; do\n' % Ndir
        #     msg += 'if [ $i -eq 0 ]; then\n'
        #     msg += '    cp $RUNDIR/$f/ATAmatrix.datATA.dat\n'
        #     msg += '    cp $RUNDIR/$f/ATd.dat .\n'
        #     msg += '    cp $RUNDIR/$f/dTd.dat .\n'
        #     msg += 'else\n'
        #     msg += '    echo $RUNDIR/$f/ATAmatrix.dat >> ATA.list\n'
        #     msg += '    echo $RUNDIR/$f/ATd.dat >> ATA.list\n'
        #     msg += '    echo $RUNDIR/$f/dTd.dat >> ATA.list\n'
        #     msg += 'fi\n'
        #     msg += 'i=$(( i+1 ))\n'
        #     msg += 'done\n'
        #     msg += '{addATA} < ATA.list\n'
        msg += 'echo "    -Adding ATA" >> progress.out\n'
        msg += 'python {addATA} $RUNDIR %s\n' % Ndir
        msg += 'cp ATA.dat ATAmatrix.dat\n'

        msg = msg.format(addATA='merge.py')

    elif setup.intype == 'FSI':
        splittingdir = setup.splittingdir
        if splittingdir.startswith('//nfs'):
            splittingdir = splittingdir.replace('//nfs/stig', '/home')
        if isinstance(setup.events, AttribDict):
            tot_ev = 0
            for segdir_name, allevents in setup.events.items():
                if tot_ev == 0:
                    events = sorted(allevents)
                    dir = '%s/%s/%s' % (splittingdir, segdir_name, events[0])
                    msg += 'echo "    -addATA" >> $RUNDIR/progress.out\n'
                    msg += 'rm ATAmatrix.dat ATd.dat dTd.dat ATA.dat\n'
                    msg += 'cp %s/ATAmatrix.dat %s/ATA.dat\n' % (dir, rundir)
                    msg += 'cp %s/ATd.dat %s/.\n' % (dir, rundir)
                    msg += 'cp %s/dTd.dat %s/.\n\n' % (dir, rundir)
                tot_ev += len(allevents)

            tot_segdir = 0
            for segdir_name, allevents in setup.events.items():
                events = sorted(allevents)
                if len(events) > 1 and tot_segdir == 0:
                    msg += '{addATA} <<EOF\n'
                    msg += '{msize}\n'
                    msg += '%d\n' % (tot_ev - 1)
                    for e in events[1:]:
                        dir = '%s/%s/%s' % (splittingdir, segdir_name, e)
                        msg += '%s/ATAmatrix.dat\n' % (dir)
                        msg += '%s/ATd.dat\n' % (dir)
                        msg += '%s/dTd.dat\n' % (dir)

                elif len(events) > 1 and tot_segdir > 0:
                    for e in events:
                        dir = '%s/%s/%s' % (splittingdir, segdir_name, e)
                        msg += '%s/ATAmatrix.dat\n' % (dir)
                        msg += '%s/ATd.dat\n' % (dir)
                        msg += '%s/dTd.dat\n' % (dir)
                tot_segdir += 1

            if len(events) > 1:
                msg += 'EOF\n'

            msg += 'cp ATA.dat ATAmatrix.dat\n'

            msg = msg.format(addATA=join(setup.bindir, setup.bins.addATA),
                             msize=setup.model_size)
    return msg


def dorun_invATA(setup):
    msg = '#-------------------#\n'
    msg += '#  Build new Model  #\n'
    msg += '#-------------------#\n'
    msg += 'echo "    -invATA" >> progress.out\n'
    msg += '{invATA} << EOF\n'
    msg += '$damp\n'
    msg += '0\n'
    msg += 'EOF\n'

    if setup.intype == 'FSI':
        # transform mnew.dat to spline format
        msg += '%s <<EOF\n' % (join(setup.bindir, setup.bins.m2sph))
        msg += '%s\n' % (setup.model.split('/')[-1])
        msg += 'EOF\n'
        msg += 'cp %s mnew_it"$iterno"_d"$damp".sph\n'
        msg = msg % (setup.model.split('/')[-1])

    msg = msg.format(invATA=join(setup.bindir, setup.bins.invATA))
    return msg


def dorun_mdcpl(setup):
    msg = '\n'
    msg += '#--------------------------#\n'
    msg += '#  Build matrix with mnew  #\n'
    msg += '#--------------------------#\n'
    if setup.intype == 'CST':
        msg += 'cp mnew.dat mcst.dat\n'
        if len(setup.modes_sc_dst) > 0:
            msg += 'echo "    -write_dst" >> progress.out\n'
            msg += '{write_dst} << EOF\n'
            for mode, smax in setup.modes_sc.items():
                msg += '%s\n' % smax
            for mode, smax in setup.modes_sc_dst.items():
                msg += '%s\n' % smax
            msg += 'EOF\n\n'
            msg = msg.format(write_dst=join(setup.bindir,
                                            setup.bins.write_dst))

        msg += 'echo "    -matrixcst" >> progress.out\n'
        msg += '{mdcpl} << EOF\n'
        for mode, smax in setup.modes_sc.items():
            # /home/deuss/bin/matrixcst
            # maximum structural deg for every mode
            msg += '0 %s\n' % smax
        if hasattr(setup, 'modes_cc'):
            # /home/deuss/bin/matrixcstC
            for mode, s in setup.modes_cc.items():
                msg += '%s %s\n' % (s[0], s[1])
        msg += 'EOF\n'
        msg = msg.format(mdcpl=join(setup.bindir, setup.bins.mdcpl))

        if len(setup.modes_sc_dst) > 0:
            msg += 'mv matrix.dat matrix-cst.dat\n\n'
            # /home/deuss/bin/matrixdst
            msg += 'echo "    -matrixdst" >> progress.out\n'
            msg += '{mdcplmq} << EOF\n'
            for mode, smax in setup.modes_sc_dst.items():
                # /home/deuss/bin/matrixcst
                # maximum structural deg for every mode
                msg += '2 %s\n' % smax
            msg += 'EOF\n'
            msg += 'mv matrix.dat matrix-dst.dat\n\n'

            msg += 'echo "    -addmdcplmq" >> progress.out\n'
            msg += '{addmdcplmq} << EOF\n'
            msg += 'matrix-cst.dat\n'
            msg += 'matrix-dst.dat\n'
            msg += 'EOF\n'
            msg = msg.format(mdcplmq=join(setup.bindir, setup.bins.mdcplmq),
                             addmdcplmq=join(setup.bindir,
                             setup.bins.addmdcplmq))

        msg += '\n'
        msg += 'echo "    -matdiag" >> progress.out\n'
        msg += '{matdiag} << EOF\n'  # /home/deuss/bin/matdiag
        msg += '0\n'
        msg += 'EOF\n\n'

        msg = msg.format(matdiag=join(setup.bindir, setup.bins.matdiag))

    if setup.intype == 'FSI':
        pbsrundir = pbsrundir_sanity_check(setup.rundir)

        msg += 'rm misfit.dat\n'
        msg += 'rm dsyni.dat dsynr.dat dobsi.dat dobsr.dat misfit-mode.dat\n'
        msg += '{mdcpl} -model {model} <<EOF\n'
        msg += '0 {s_deg}\n'
        msg += '{rotation}\n'
        msg += '{ellip}\n'
        msg += 'EOF\n\n'
        msg += '{matdiag} << EOF\n'  # /home/deuss/bin/matdiag
        msg += '0\n'
        msg += 'EOF\n\n'
        msg = msg.format(mdcpl=join(setup.bindir, setup.bins.mdcpl),
                         matdiag=join(setup.bindir, setup.bins.matdiag),
                         s_deg=setup.s_deg,
                         model=join(pbsrundir, setup.model.split('/')[-1]),
                         rotation=setup.rotation,
                         ellip=setup.ellip)
    return msg


def dorun_synseis_inv(setup):
    msg = '#----------------------------#\n'
    msg += '#  Calculate new synthetics  #\n'
    msg += '#----------------------------#\n'
    msg += 'rm misfit.dat\n'
    msg += 'rm dsyni.dat dsynr.dat dobsi.dat dobsr.dat misfit-mode.dat\n'

    if setup.intype == 'CST':
        msg += 'echo "    -synseis_inv" >> progress.out\n'
    elif setup.intype == 'FSI':
        msg += get_segments(list(setup.seg_suffix.values())[0])
    msg += 'for FILE in "${{SEGMENTS[@]}}"\n'
    msg += 'do\n'
    msg += '{synseis_inv} '
    msg += get_FILE_extension(list(setup.seg_suffix.values())[0])
    msg += ' << EOF\n'
    msg += 'rr.dat\n'
    msg += 'rr1.dat\n'
    msg += 'ww.dat\n'
    msg += 'EOF\n'
    msg += 'done\n\n'
    msg = msg.format(synseis_inv=join(setup.bindir, setup.bins.synseis_inv))

    msg += 'count=`grep "NaN" "$EVENT".misfit | wc -l`\n'
    msg += '\n'
    msg += 'if [ $count != 0 ]; then\n'
    msg += '        stat=`grep -B 3 "NaN" "$EVENT".misfit | '
    msg += "cut -d'-' -f1 | cut -d' ' -f1`\n"
    msg += '        echo "synseis    -NaN found in Station: '
    msg += '$SEGDIR/$EVENT $stat" >> $MAINDIR/progress.err\n'
    msg += '        sed -i "/$damp $SEGDIR $EVENT/s/$/ NaN found in Station: '
    msg += '$SEGDIR/$EVENT $stat/" $MAINDIR/progress.out\n'
    msg += '        #rm misfit.dat\n'
    msg += '        #rm dsyni.dat dsynr.dat dobsi.dat '
    msg += 'dobsr.dat misfit-mode.dat\n'
    msg += '        cd $MAINDIR\n'
    msg += '        exit 404\n'
    msg += 'fi\n'
    return msg


def dorun_synseis_inv_array(setup, uncertainties=False):
    msg = '#------------------------#\n'
    msg += '#  Calculate synthetics  #\n'
    msg += '#-----------------------#\n'
    msg += 'rm misfit.dat 2> /dev/null\n'
    msg += 'rm dsyni.dat dsynr.dat dobsi.dat dobsr.dat '
    msg += 'misfit-mode.dat 2> /dev/null\n\n'

    if setup.intype == 'CST':
        if len(setup.events) > 1:
            msg += 'PO=$MAINDIR/progress.out\n'
            msg += 'sed -i "/$iterno $damp $SEGDIR $EVENT/s/$/ '
            msg += 'synseis_inv/" $PO\n'
        else:
            msg += 'PO=$MAINDIR/progress.out\n'
            msg += 'sed -i "/$iterno $damp $EVENT/s/$/ synseis_inv/" $PO\n'

    msg += '{synseis_inv} $EVENT << EOF\n'
    msg += 'rr.dat\n'
    msg += 'rr1.dat\n'
    msg += 'ww.dat\n'
    msg += 'EOF\n\n'
    if uncertainties is False:
        msg += synseis_sanity_check(setup)

    msg = msg.format(synseis_inv=join(setup.bindir, setup.bins.synseis_inv))
    return msg


def synseis_sanity_check(setup):
    msg = 'count=`grep "NaN" "$EVENT".misfit | wc -l`\n\n'
    msg += 'if [ $count != 0 ]; then\n'
    msg += '\tstat=`grep -B 3 "NaN" "$EVENT".misfit'
    msg += " | cut -d'-' -f1 | cut -d' ' -f1`\n"
    msg += '\techo "synseis    -NaN found in Station: $SEGDIR/$EVENT '
    msg += '$stat" >> $MAINDIR/progress.err\n'
    if len(setup.events) > 1:
        msg += '\tsed -i "/$iterno $damp $SEGDIR $EVENT/s/$/ '
        msg += 'NaN found in Station: '
        msg += '$SEGDIR/$EVENT $stat/" $PO\n'
    else:
        msg += '\tsed -i "/$iterno $damp $EVENT/s/$/ NaN found in Station: '
        msg += '$EVENT $stat/" $PO\n'
    msg += '\trm misfit.dat 2> /dev/null\n'
    msg += '\trm dsyni.dat dsynr.dat dobsi.dat dobsr.dat '
    msg += 'misfit-mode.dat 2> /dev/null\n'
    msg += '\tcd $MAINDIR\n'
    msg += '\texit 404\n'
    msg += 'fi\n\n'

    msg = msg.format(synseis_inv=join(setup.bindir, setup.bins.synseis_inv))
    return msg


def dorun_misfits(setup):
    msg = '\n'
    msg += '#--------------------#\n'
    msg += '#  Calculate misfit  #\n'
    msg += '#--------------------#\n'
    if setup.intype == 'FSI':
        msg += '{avmisfit_allmodes}\n'
    msg += 'echo "    -avmisfit" >> progress.out\n'
    msg += '{avmisfit}\n'  # /home/deuss/bin/avmisfit
    msg += 'echo "    -posteriorCd" >> progress.out\n'
    msg += '{posterioriCd}\n'  # /home/deuss/bin/posterioriCd
    msg += 'echo "    -inversion_summary" >> progress.out\n'
    msg += 'python {invout} $RUNDIR $iterno {summary_dir} {schema}\n'
    msg += 'rm inversion.out 2> /dev/null\n'
    # #### dorun-newC3 new end
    if setup.intype == 'CST':
        msg += 'cp mnew-error.dat mnew.dat\n'
        msg += 'cp mnew.dat {summary_dir}/mnew-it"$iterno"-d"$damp".dat\n'
        msg += 'cp mnew.dat mold.dat\n'
        msg += 'cp mnew.dat mcst.dat\n'
    elif setup.intype == 'FSI':
        msg += 'cp misfit.dat misfit_it${{iterno}}_d${{damp}}.dat\n'
        msg += 'cp mnew.dat mnew_it${{iterno}}_d${{damp}}.dat\n'
        msg += 'cp mnew.dat mold.dat\n'

    msg = msg.format(avmisfit_allmodes=join(setup.bindir,
                                            setup.bins.avmisfit_allmodes),
                     avmisfit=join(setup.bindir, setup.bins.avmisfit),
                     posterioriCd=join(setup.bindir, setup.bins.posterioriCd),
                     invout="write_inversion_out.py",
                     summary_dir=setup.inversion_outdir,
                     schema=setup.schema)
    return msg


def dorun_misfits_array(setup, uncertainties=False):
    if uncertainties is True:
        setup.inversion_outdir = 'uncertainties'
    # In this case we use the same function for the startup and new runs
    msg = '\n'
    msg += '#--------------------#\n'
    msg += '#  Calculate misfit  #\n'
    msg += '#--------------------#\n'
    msg += 'echo "    -avmisfit" >> progress.out\n'
    if len(setup.events) > 1:
        msg += '{avmisfit_allmodes}\n'
    msg += '{avmisfit}\n'  # /home/deuss/bin/avmisfit
    if len(setup.segmentsdir) > 1:
        for idir, s in enumerate(setup.segmentsdir.keys()):
            msg += 'python {invout} $RUNDIR $(($iterno-1)) '
            if uncertainties is True:
                msg2 = '{summary_dir} "%s-$set_filename"' % s
            else:
                msg2 = '{summary_dir} %s' % s
            if idir != 0:
                msg2 += ' no_summary'
            msg += msg2
            msg += '\n'
    else:
        msg += 'python {invout} $RUNDIR $(($iterno-1)) {summary_dir}\n'

    msg += '\n'

    msg += 'if [[ ! $iterno -eq 1 ]]; then\n'
    msg += 'echo "    -posteriorCd" >> progress.out\n'
    msg += '{posterioriCd}\n'  # /home/deuss/bin/posterioriCd
    msg += 'echo "    -inversion_summary" >> progress.out\n'
    msg += 'rm inversion.out 2> /dev/null\n'
    msg += 'cp mnew-error.dat mnew.dat\n'
    if setup.intype == 'CST':
        msg += 'cp mnew.dat mcst.dat\n\n'

    if uncertainties is False:
        msg += 'cp mnew.dat {summary_dir}/mnew-it"$(($iterno-1))"'
        msg += '-d"$damp".dat\n'
    msg += 'cp mnew.dat mold.dat\n\n'
    msg = msg.format(avmisfit_allmodes=join(setup.bindir,
                                            setup.bins.avmisfit_allmodes),
                     avmisfit=join(setup.bindir, setup.bins.avmisfit),
                     posterioriCd=join(setup.bindir, setup.bins.posterioriCd),
                     invout="write_inversion_out.py",
                     summary_dir=setup.inversion_outdir)
    return msg


def dorun_summary(setup, uncertainties=False):
    msg = '\n'
    msg += '#---------------------#\n'
    msg += '#  INVERSION SUMMARY  #\n'
    msg += '#---------------------#\n'
    if uncertainties is True:
        msg += 'cp mnew.dat {summary_dir}/mnew-it"$(($iterno-1))"'
        msg += '-d"$damp"-"$set_filename".dat\n'
    else:
        msg += 'cp Cmatrix.dat {summary_dir}/'
        msg += 'Cmatrix-it"$(($iterno-1))"-d"$damp".dat\n'
        msg += 'cp Rmatrix.dat {summary_dir}/'
        msg += 'Rmatrix-it"$(($iterno-1))"-d"$damp".dat\n'
    if setup.keep_all_output is not True:
        msg += 'for mf in *.misfit; do cp $mf '
        msg += '{summary_dir}/${{mf%.misfit}}_d"$damp".misfit; done\n'
    msg = msg.format(summary_dir=setup.inversion_outdir)
    return msg


def dorun_cleanup(setup):
    msg = 'rm inversion.dat *misfit* mdcpl.out *.bin mnew* Cmat* Rmat*\n'
    if setup.keep_all_output is not True:
        msg += 'ls *dat |  xargs -I {{}} find {{}} -type f -not -name '
        msg += "'???????.dat' -delete \n"
    return msg


def dorun_SLURM_header(setup):
    pbsrundir = pbsrundir_sanity_check(setup.rundir)
    msg = '#!/bin/bash\n'
    msg += '#SBATCH --exclusive\n'
    msg += '#SBATCH --mem=200G\n'
    msg += '#SBATCH --partition=%s\n' % setup.partition
    msg += '#SBATCH --nodes=1\n'
    msg += '#SBATCH --cpus-per-task=48\n'
    msg += '#SBATCH --threads-per-core=2\n'
    if setup.walltime is not None:
        msg += '#SBATCH --time=%s\n' % setup.walltime
    else:
        msg += '#SBATCH --time=24:00:00\n'
    msg += '#SBATCH --job-name=%s.job\n' % setup.pbsjob
    # msg += '#SBATCH --output=%s/main.%%A.%%a.out\n' % (pbsrundir)
    msg += '#SBATCH --error=%s/main.%%A.%%a.err\n' % (pbsrundir)
    msg += '#SBATCH --export=all\n'

    if hasattr(setup, 'pbsmail'):
        msg += '#SBATCH --mail-type=ALL\n'
        msg += '#SBATCH --mail-user=%s\n' % setup.pbsmail
    msg += 'RUNDIR=%s\n' % pbsrundir
    return msg


def dorun_submission_header(setup, type='default', uncertainties=False):

    pbsrundir = pbsrundir_sanity_check(setup.rundir)

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        msg = '#!/bin/bash\n'

        if type in ('default', 'submission', 'depend'):
            msg += '#SBATCH --exclusive\n'
            msg += '#SBATCH --mem=200G\n'

            if setup.exclude is not None:
                msg += '#SBATCH --exclude=%s\n' % setup.exclude

            msg += '#SBATCH --partition=%s\n' % setup.partition

            if hasattr(setup, 'pbsmail'):
                msg += '#SBATCH --mail-type=ALL\n'
                msg += '#SBATCH --mail-user=%s\n' % setup.pbsmail
            msg += '#SBATCH --nodes=%s\n' % setup.nodes
            msg += '#SBATCH --cpus-per-task=48\n'
            msg += '#SBATCH --threads-per-core=2\n'
            if setup.walltime is not None:
                msg += '#SBATCH --time=%s\n' % setup.walltime
            else:
                msg += '#SBATCH --time=24:00:00\n'

        if type == 'default':
            if uncertainties is False:
                msg += '#SBATCH --job-name=%s\n' % setup.pbsjob
                msg += '#SBATCH --error=%s/main.%%A.%%a.err\n' % (pbsrundir)
            else:
                msg += '#SBATCH --job-name=%s.uncertainties\n' % setup.pbsjob
                msg += '#SBATCH --error=%s/uncert.%%A.%%a.err\n' % (pbsrundir)
            msg += '#SBATCH --export=all\n'
            if setup.rundir.startswith('/quanta'):
                msg += 'MOM=%s\n' % pbsrundir
                msg += 'RUNDIR=/tmp/$USER\n'
                msg += 'if [ -e $RUNDIR ]; then rm -rf $RUNDIR; fi\n'
                msg += 'mkdir -p $RUNDIR\n'
                msg += 'rsync -r $MOM/* $RUNDIR/.\n'
            else:
                msg += 'RUNDIR=%s\n' % pbsrundir

        elif type == 'submission':
            msg += '#SBATCH --job-name=%s.job\n' % setup.pbsjob
            # msg += '#SBATCH --output=%s/log/main.%%A.%%a.out\n' % (pbsrundir)
            msg += '#SBATCH --error=%s/log/main.%%A.%%a.err\n' % (pbsrundir)
            msg += '#SBATCH --export=all\n'
            msg += 'RUNDIR=%s\n' % pbsrundir

        elif type == 'depend':
            # msg += '#SBATCH --output=%s/log/array.%%A.%%a.out\n' % pbsrundir
            msg += '#SBATCH --error=%s/log/array.%%A.%%a.err\n' % pbsrundir
            msg += 'RUNDIR=%s\n' % pbsrundir

        elif type == 'array':
            if not setup.rundir.startswith('/quanta'):
                msg += 'MAINDIR=%s\n' % pbsrundir
            if len(setup.events) > 1:
                msg += 'EVENT=$1\n'
                msg += 'iterno=$2\n'
                msg += 'ed_it=$3\n'
                msg += 'SEGDIR=$4\n'
                msg += 'damp=$5\n'
                msg += 'MAINDIR=$6\n\n'
                if setup.intype == 'CST':
                    msg += 'RUNDIR=$MAINDIR/$SEGDIR/$EVENT\n'
            else:
                msg += 'EVENT=$1\n'
                msg += 'iterno=$2\n'
                msg += 'ed_it=$3\n'
                msg += 'damp=$4\n'
                msg += 'MAINDIR=$5\n\n'
                if setup.intype == 'CST':
                            msg += 'RUNDIR=$MAINDIR/$EVENT\n'

            if setup.intype == 'FSI':
                splittingdir = setup.splittingdir
                if splittingdir.startswith('//nfs'):
                    splittingdir = splittingdir.replace('//nfs/stig/',
                                                        '/home/')
                msg += 'SPLITDIR=%s\n' % (splittingdir)

                if len(setup.events) > 1:
                    msg += 'RUNDIR=$SPLITDIR/$SEGDIR/$EVENT\n'

                else:
                    msg += 'RUNDIR=$SPLITDIR/$EVENT\n'

    else:
        msg = ''
        if hasattr(setup, 'pbsmail'):
            msg += '#PBS -m abe -M %s\n' % setup.pbsmail
        if type == 'default':
            msg += '#PBS -N %s\n' % setup.pbsjob
            msg += '#PBS -l nodes=1:ppn=1\n'
            msg += '#PBS -e log.$PBS_JOBNAME.$PBS_JOBID.err\n'
            msg += '#PBS -o log.$PBS_JOBNAME.$PBS_JOBID.out\n'
            msg += 'RUNDIR=%s\n' % pbsrundir

        elif type == 'depend':
            msg += '#PBS -l nodes=1:ppn=1\n'
            msg += '#PBS -e log/log.$PBS_JOBNAME.$PBS_JOBID.err\n'
            msg += '#PBS -o log/log.$PBS_JOBNAME.$PBS_JOBID.out\n'
            msg += 'RUNDIR=%s\n' % pbsrundir

        elif type == 'array':
            if (
                isinstance(setup.segmentsdir, AttribDict) or
                isinstance(setup.segmentsdir, dict)
               ):
                tot_ev = 0
                for segdir, allevents in setup.events.items():
                    tot_ev += len(allevents)
                n = tot_ev

                if n > 96:
                    t = 96
                else:
                    t = n
            else:
                n = len(setup.events)
                t = n
            msg += '#PBS -l walltime=1:00:00,nodes=1:ppn=1\n'
            msg += '#PBS -t 1-%s%%%s\n' % (n, t)
            msg += '#PBS -e log/log.$PBS_JOBNAME.$PBS_JOBID.err\n'
            msg += '#PBS -o log/log.$PBS_JOBNAME.$PBS_JOBID.out\n\n'
            msg += 'MAINDIR=%s\n' % pbsrundir

    return msg


def script(setup, output='pbs'):
    """
    Creates a pbs or bash script to submit an inversion, using the information
    given in the setup object.
    """

    if setup.intype == 'CST':
        if output in [None, 'pbs', 'SLURM', 'slurm']:
            msg = dorun_submission_header(setup)
        elif output == 'bash':
            msg = 'RUNDIR=%s\n' % setup.rundir

        msg += 'cd $RUNDIR\n'

        msg += '\n'
        msg += '#----------------#\n'
        msg += '#  DAMPING LOOP  #\n'
        msg += '#----------------#\n'
        msg += 'for damp in %s\n' % get_damping_string(setup)
        msg += 'do\n'
        msg += 'echo "Starting run for damping: $damp" >> progress.out\n'
        # /home/deuss/modes/allsplitting/dorun-itmC3.sh ${lmax} 4 0 ${damp}

        msg += dorun_startup(setup)

        msg += '\n'
        msg += '#------------------------#\n'
        msg += '#  ITERATION LOOP START  #\n'
        msg += '#------------------------#\n'
        # ##/home/deuss/modes/allsplitting/dorun-derivC3.sh $1 $2 $3\n'
        m = 'for (( iterno=%s; iterno<=%s; iterno++ ))\n'
        msg += m % (setup.iterations[0], setup.iterations[1])
        msg += 'do\n'
        msg += 'echo "  Iteration: $iterno" >> progress.out\n'

        msg += dorun_cst_partials(setup)

        if len(setup.modes_sc_dst) > 0:
            msg += dorun_dst_partials(setup)

        msg += dorun_buildATA(setup)

        msg += dorun_invATA(setup)

        # '/home/deuss/modes/allsplitting/dorun-newC3.sh $1 $2 $3\n'
        msg += dorun_mdcpl(setup)

        msg += dorun_synseis_inv(setup)

        msg += dorun_misfits(setup)

        msg += 'done\n'
        msg += '\n'
        msg += '#----------------------#\n'
        msg += '#  ITERATION LOOP END  #\n'
        msg += '#----------------------#\n'

        msg += dorun_summary(setup)

        msg += 'done\n'
        msg += '#--------------------#\n'
        msg += '#  DAMPING LOOP END  #\n'
        msg += '#--------------------#\n'
        msg += dorun_cleanup(setup)
        msg += 'echo "Finished Run" >> progress.out\n'

    elif setup.intype == 'FSI':
        msg = dorun_submission_header(setup)
        msg += 'cd $RUNDIR\n'

        modes = ''
        for mode in setup.modedirs:
            mdir = join(setup.splittingdir, mode)
            modes = modes + "%s " % (mdir)

        msg += dorun_preprundir(setup)

        msg += '\n'
        msg += '#----------------#\n'
        msg += '#  DAMPING LOOP  #\n'
        msg += '#----------------#\n'
        msg += 'for damp in %s\n' % get_damping_string(setup)
        msg += 'do\n'
        msg += 'echo "Starting run for damping: $damp" >> progress.out\n'

        msg += 'echo "Starting up" >> progress.out\n'
        msg += 'for mdir in %s\n' % (modes)
        msg += 'do\n'
        msg += 'cd $mdir\n'
        msg += dorun_startup(setup)
        msg += 'done\n'

        msg += '\n'
        msg += '#------------------------#\n'
        msg += '#  ITERATION LOOP START  #\n'
        msg += '#------------------------#\n'
        msg += 'for (( iterno=%s; iterno<=%s; iterno++ ))\n'
        msg = msg % (setup.iterations[0], setup.iterations[1])
        msg += 'do\n'
        msg += 'echo "  Iteration: $iterno" >> progress.out\n'

        msg += 'echo "    -fsi partials" >> progress.out\n'
        msg += 'echo "    -buildATA" >> progress.out\n'
        msg += '# BEGIN LOOP THROUGH MODE FOLDERS\n'
        msg += 'for mdir in %s\n' % (modes)
        msg += 'do\n'
        msg += 'cd $mdir\n'
        msg += dorun_fsi_partials(setup)

        msg += dorun_buildATA(setup)

        msg += 'done\n'  # end loop over mode folders
        msg += '# END LOOP THROUGH MODE FOLDERS\n'

        # GO BACK TO INVERSION FOLDER
        msg += '#\n'
        msg += 'cd $RUNDIR\n'

        msg += dorun_addATA(setup)

        msg += dorun_invATA(setup)

        # GO TO EACH MODES-FOLDER TO CALCULATE MDCPL, MATDIAG, SYNSEIS, MISFIT
        msg += '#\n'
        msg += '# BEGIN LOOP THROUGH MODE FOLDERS\n'
        msg += 'for mdir in %s\n' % (modes)
        msg += 'do\n'
        msg += 'cd $mdir\n'

        msg += dorun_mdcpl(setup)

        msg += dorun_synseis_inv(setup)

        msg += 'done\n'  # end loop over mode folders
        msg += '# END LOOP THROUGH MODE FOLDERS\n'
        msg += '#\n'

        # AND BACK TO INVERSION FOLDER
        msg += 'cd $RUNDIR\n'

        msg += dorun_misfits(setup)

        msg += 'done\n'
        msg += '\n'
        msg += '#----------------------#\n'
        msg += '#  ITERATION LOOP END  #\n'
        msg += '#----------------------#\n'

        msg += dorun_summary(setup)

        msg += 'done\n'
        msg += '#--------------------#\n'
        msg += '#  DAMPING LOOP END  #\n'
        msg += '#--------------------#\n'
        msg += dorun_cleanup(setup)
        msg += 'echo "Finished Run" >> progress.out\n'

    if output == 'bash':
        pbs_file = 'run_me.sh'
    elif output == 'slurm':
        pbs_file = 'sbatch_me.sh'
    else:
        pbs_file = 'qsub_me.pbs'

    if setup.rundir.startswith('/scratch'):
        src = setup.rundir
        dst = setup.rundir.replace('/scratch', '/quanta1/home')
        dst = '/'.join(dst.split('/')[:-1])
        msg += 'SRC=%s\n' % (src)
        msg += 'DST=%s\n' % (dst)
        msg += 'mkdir -p %s \n' % dst
        msg += 'rsync -r --delete $SRC $DST'

    pbs_path = join(setup.rundir, pbs_file)
    with open(pbs_path, 'w+') as fh:
        fh.write(msg)
    return


def startup_script(setup):
    # for array pbs dependencies
    """
    Creates a pbs script to submit an inversion, using the information
    given in the setup object.
    """
    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        if setup.intype == 'CST':
            msg = '#!/bin/bash\n'
            msg += 'damp=${{1}}\n'
            msg += 'RUNDIR=${{2}}\n\n'
        elif setup.intype == 'FSI':
            msg = '#!/bin/bash\nsdir=${{1}}\n'
            msg += 'RUNDIR=%s\n' % setup.rundir
    else:
        msg = dorun_submission_header(setup, type='depend')

    if setup.intype == 'CST':
        m = 'cd $RUNDIR\n'
        m += '\n'

        # Startmodel mzero.dat and mcst.dat are generated by create.startmodel
        m += '#---------------------------#\n'
        m += '#  STARTUP WITH ZERO MODEL  #\n'
        m += '#---------------------------#\n'

        m += 'echo "Starting run for damping: $damp" >> progress.out\n'
        m += 'echo "Starting up" >> progress.out\n'
        m += 'rm inversion.dat inversion.out att.dat 2> /dev/null\n'
        m += '\n'
        m += 'cp mzero.dat mcst.dat\n'
        m += 'cp mzero.dat inversion_out/.\n\n'
        m += 'echo "    -rot+ellip" >> progress.out\n'
        m += '{rot_ellip}\n\n'
        m = m.format(rot_ellip=join(setup.bindir, setup.bins.rot_ellip))

        if len(setup.modes_sc_dst) > 0:
            m += 'echo "    -write_dst" >> progress.out\n'
            m += '{write_dst} << EOF\n'
            for mode, smax in setup.modes_sc.items():
                m += '%s\n' % smax
            for mode, smax in setup.modes_sc_dst.items():
                m += '%s\n' % smax
            m += 'EOF\n\n'
            m = m.format(write_dst=join(setup.bindir, setup.bins.write_dst))

        m += 'echo "    -matrixcst" >> progress.out\n'
        m += '{mdcpl} << EOF\n'
        for mode, smax in setup.modes_sc.items():
            # /home/deuss/bin/matrixcst
            # maximum structural deg for every mode
            m += '0 %s\n' % smax
        if hasattr(setup, 'modes_cc'):
            # /home/deuss/bin/matrixcstC
            for mode, s in setup.modes_cc.items():
                m += '%s %s\n' % (s[0], s[1])
        m += 'EOF\n'
        m = m.format(mdcpl=join(setup.bindir, setup.bins.mdcpl))

        if len(setup.modes_sc_dst) > 0:
            m += 'mv matrix.dat matrix-cst.dat\n\n'
            # /home/deuss/bin/matrixdst
            m += 'echo "    -matrixdst" >> progress.out\n'
            m += '{mdcplmq} << EOF\n'
            for mode, smax in setup.modes_sc_dst.items():
                # /home/deuss/bin/matrixcst
                # maximum structural deg for every mode
                m += '2 %s\n' % smax
            m += 'EOF\n'
            m += 'mv matrix.dat matrix-dst.dat\n\n'
            m += 'echo "    -addmdcplmq" >> progress.out\n'
            m += '{addmdcplmq} << EOF\n'
            m += 'matrix-cst.dat\n'
            m += 'matrix-dst.dat\n'
            m += 'EOF\n\n'
            m = m.format(mdcplmq=join(setup.bindir, setup.bins.mdcplmq),
                         addmdcplmq=join(setup.bindir, setup.bins.addmdcplmq))

    elif setup.intype == 'FSI':
        m = 'cd $sdir\n'
        m += 'echo "    -mdcpl" >> $RUNDIR/progress.out\n'
        m += '{mdcpl} -model {zeromodel} <<EOF\n'
        m += '0 {s_deg}\n'
        m += '{rotation}\n'
        m += '{ellip}\n'
        m += 'EOF\n\n'
        m = m.format(mdcpl=join(setup.bindir, setup.bins.mdcpl),
                     zeromodel=setup.model, rotation=setup.rotation,
                     ellip=setup.ellip, s_deg=setup.s_deg)

    msg += m
    msg += '\n'
    msg += 'echo "    -matdiag" >> $RUNDIR/progress.out\n'
    msg += '{matdiag} << EOF\n'
    msg += '0\n'
    msg += 'EOF\n\n'
    msg = msg.format(matdiag=join(setup.bindir, setup.bins.matdiag))

    if setup.intype == 'CST':
        msg += 'cp mcst.dat mold.dat\n\n'
    elif setup.intype == 'FSI':
        pass

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        pbs_file = 'sbatch_startup.sh'
    else:
        pbs_file = 'qsub_startup.pbs'

    pbs_path = join(setup.rundir, pbs_file)
    with open(pbs_path, 'w+') as fh:
        fh.write(msg)
    return


def array_script(setup, uncertainties=False):
    """
    Creates a pbs script to submit an inversion, using the information
    given in the setup object.
    """

    msg = dorun_submission_header(setup, type='array', uncertainties=False)

    if setup.intype == 'CST':
        msg += 'if [ ! -d $RUNDIR ]; then\n'
        msg += '        echo "$RUNDIR does not exist, exiting array"\n'
        msg += '        exit\n'
        msg += 'fi\n'
        msg += 'cp $MAINDIR/rr*.dat $RUNDIR\n'
        msg += 'cp $MAINDIR/ww.dat $RUNDIR\n'
        msg += 'cd $RUNDIR\n'
        msg += '\n'
    elif setup.intype == 'FSI':
        msg += 'cp $SPLITDIR/$SEGDIR/rr*.dat $RUNDIR\n'
        msg += 'cp $SPLITDIR/$SEGDIR/ww.dat $RUNDIR\n'
        msg += 'cd $RUNDIR\n'
        msg += '\n'

    msg += dorun_synseis_inv_array(setup, uncertainties)

    # No partials in last iteration
    msg += 'if [[ $iterno -ne $ed_it ]]; then\n'
    if setup.intype == 'CST':
        msg += dorun_cst_partials_array(setup)
        if len(setup.modes_sc_dst) > 0:
            msg += dorun_dst_partials_array(setup)
        msg += dorun_buildATA_array(setup)
    elif setup.intype == 'FSI':
        msg += dorun_fsi_partials_array(setup)
    msg += '\n'
    msg += 'fi\n\n'

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):

        if uncertainties is False:
            pbs_file = 'sbatch_array.sh'
        else:
            pbs_file = 'uncertainties/sbatch_array.sh'
    else:
        pbs_file = 'qsub_array.pbs'

    pbs_path = join(setup.rundir, pbs_file)
    with open(pbs_path, 'w+') as fh:
        fh.write(msg)
    os.chmod(pbs_path, 0o777)
    return


def new_script(setup, uncertainties=False):
    """
    Creates a pbs script to submit an inversion, using the information
    given in the setup object.
    """
    def cat_misfits(setup, uncertainties):
        msg = "while read event; do\n"
        msg += "    weight=`echo $event | awk '{print $NF}'`\n"
        msg += "    if [[ $weight =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then\n"
        msg += "        EVENT=`echo $event | awk '{print $2}'`\n"
        msg += "    else\n"
        msg += "        EVENT=`echo $event | awk '{print $NF}'`\n"
        msg += "    fi\n"
        msg += "    SEGDIR=`echo $event | awk '{print $1}'`\n"
        if len(setup.events) > 1:
            msg += '\tif [ -e $SEGDIR/$EVENT/$EVENT.misfit ]; then\n'
            msg += '\t\tcat $SEGDIR/$EVENT/dobsr.dat >> dobsr.dat\n'
            msg += '\t\tcat $SEGDIR/$EVENT/dobsi.dat >> dobsi.dat\n'
            msg += '\t\tcat $SEGDIR/$EVENT/dsynr.dat >> dsynr.dat\n'
            msg += '\t\tcat $SEGDIR/$EVENT/dsyni.dat >> dsyni.dat\n'
            msg += '\t\tcat $SEGDIR/$EVENT/misfit.dat >> misfit-$SEGDIR.dat\n'
            msg += '\t\tcat $SEGDIR/$EVENT/$EVENT.misfit >> $EVENT.misfit\n'
            msg += '\telse\n'
            msg += '\t\techo "No file: $damp, $iterno, $SEGDIR/$EVENT" '
            msg += '>> $RUNDIR/progress.err\n'
            msg += '\tfi\n'
        else:
            msg += '\tif [ -e $RUNDIR/$EVENT/$EVENT.misfit ]; then\n'
            msg += '\t\tcat $RUNDIR/$EVENT/dobsr.dat >> dobsr.dat\n'
            msg += '\t\tcat $RUNDIR/$EVENT/dobsi.dat >> dobsi.dat\n'
            msg += '\t\tcat $RUNDIR/$EVENT/dsynr.dat >> dsynr.dat\n'
            msg += '\t\tcat $RUNDIR/$EVENT/dsyni.dat >> dsyni.dat\n'
            msg += '\t\tcat $RUNDIR/$EVENT/misfit.dat >> misfit.dat\n'
            msg += '\telse\n'
            msg += '\t\techo "No file: $damp, $iterno, $EVENT" '
            msg += '>> $RUNDIR/progress.err\n'
            msg += '\tfi\n'
        if uncertainties is True:
            msg += 'done < $set_name\n\n'
        else:
            msg += 'done < allevents\n\n'
        return msg

    if setup.intype == 'FSI':
        splittingdir = setup.splittingdir
        if splittingdir.startswith('//nfs'):
            splittingdir = splittingdir.replace('//nfs/stig/', '/home/')

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        msg = '#!/bin/bash\n'
        msg += 'damp=${1}\n'
        msg += 'iterno=${2}\n'
        msg += 'ed_it=${3}\n'

        if uncertainties is True:
            msg += 'set_name=${4}\n'
            msg += 'RUNDIR=${5}\n'
            msg += 'set_filename=`echo $set_name | '
            msg += "rev | cut -d'/' -f1 | rev`\n"
        else:
            msg += 'RUNDIR=${4}\n'

    else:
        msg = dorun_submission_header(setup, type='depend')

    msg += 'cd $RUNDIR\n'
    msg += 'rm dsyni.dat dsynr.dat dobsi.dat dobsr.dat 2> /dev/null\n'
    msg += 'rm misfitfiles.in misfit*dat *.misfit 2> /dev/null\n'
    msg += '\n'
    if setup.intype == 'FSI':
        msg += 'SPLITDIR=%s\n' % (splittingdir)
        msg += 'cd $SPLITDIR\n'
    msg += '#----------------------#\n'
    msg += '#  Concatenate events  #\n'
    msg += '#----------------------#\n'
    msg += 'echo "    -Concatenate events" >> $RUNDIR/progress.out\n'

    if setup.intype == 'CST':
        msg += cat_misfits(setup, uncertainties)
        msg += nan_check_misfits()
        if len(setup.events) > 1:
            for segname, sdir in setup.segmentsdir.items():
                msg += 'N=`cat misfit-%s.dat | wc -l`\n' % segname
                msg += 'if [ $N -gt 0 ]; then\n'
                msg += '\techo misfit-%s.dat >> misfitfiles.in \n' % segname
                msg += 'fi\n'

    elif setup.intype == 'FSI':
        # maybe also cat_misfits here?
        msg += 'while IFS=" " read -r SEGDIR EVENT; do\n'
        msg += 'cat $SEGDIR/$EVENT/dobsr.dat >> $RUNDIR/dobsr.dat\n'
        msg += 'cat $SEGDIR/$EVENT/dobsi.dat >> $RUNDIR/dobsi.dat\n'
        msg += 'cat $SEGDIR/$EVENT/dsynr.dat >> $RUNDIR/dsynr.dat\n'
        msg += 'cat $SEGDIR/$EVENT/dsyni.dat >> $RUNDIR/dsyni.dat\n'
        msg += 'cat $SEGDIR/$EVENT/misfit.dat >> '
        msg += '$RUNDIR/misfit-$SEGDIR.dat\n'
        msg += 'cat $SEGDIR/$EVENT/$EVENT.misfit >> '
        msg += '$RUNDIR/$EVENT.misfit\n'
        msg += 'done < $RUNDIR/allevents\n\n'
        for segname, sdir in setup.segmentsdir.items():
            msg += 'echo misfit-%s.dat >> ' % segname
            msg += '$RUNDIR/misfitfiles.in \n'
        msg += 'cd $RUNDIR\n'

    msg += dorun_misfits_array(setup, uncertainties)
    msg += 'fi\n'
    msg += '\n'
    msg += 'if [[ $iterno -ne $ed_it ]]; then\n'
    msg += 'echo "  End of iteration" $(($iterno-1)) >> progress.out\n\n'
    msg += dorun_addATA(setup, uncertainties)
    msg += dorun_invATA(setup)

    if setup.intype == 'CST':
        msg += dorun_mdcpl(setup)
    elif setup.intype == 'FSI':
        segdirs = ''
        for segname, sdir in setup.segmentsdir.items():
            if sdir.startswith('//nfs'):
                sdir = sdir.replace('//nfs/stig/', '/home/')
            segdirs = segdirs + "%s " % sdir
        msg += 'for sdir in %s\n' % segdirs
        msg += 'do\n'
        msg += 'cd $sdir\n'
        msg += dorun_mdcpl(setup)
        msg += 'done\n'

    msg += 'else\n'
    msg += dorun_summary(setup, uncertainties)

    if setup.intype == 'CST' and len(setup.segmentsdir) > 1:
        # if uncertainties is True:
        #     msg += 'while IFS=" " read -r SEGDIR EVENT WEIGHT; do\n'
        # else:
        #     msg += 'while IFS=" " read -r SEGDIR EVENT; do\n'
        msg += "while read event; do\n"
        msg += "    weight=`echo $event | awk '{print $NF}'`\n"
        msg += "    if [[ $weight =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then\n"
        msg += "        EVENT=`echo $event | awk '{print $2}'`\n"
        msg += "    else\n"
        msg += "        EVENT=`echo $event | awk '{print $NF}'`\n"
        msg += "    fi\n"
        msg += "    SEGDIR=`echo $event | awk '{print $1}'`\n"
        msg += '    cat $SEGDIR/$EVENT/$EVENT.misfit >> '
        msg += '$SEGDIR/$EVENT-d$damp.misfit\n'

        if uncertainties is True:
            msg += 'done < $set_name\n\n'
        else:
            msg += 'done < allevents\n\n'
    elif setup.intype == 'FSI' and len(setup.segmentsdir) > 1:
        msg += 'while IFS=" " read -r SEGDIR EVENT; do\n'
        msg += 'cat $SEGDIR/$EVENT/$EVENT.misfit >> '
        msg += '$SEGDIR/$EVENT-d$damp.misfit\n'
        msg += 'done < $RUNDIR/allevents\n\n'
    msg += 'echo "  End of inversion" >> progress.out\n\n'
    msg += 'fi\n\n'

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        if uncertainties is True:
            pbs_file = 'uncertainties/sbatch_new.sh'
        else:
            pbs_file = 'sbatch_new.sh'
    else:
        pbs_file = 'qsub_new.pbs'

    pbs_path = join(setup.rundir, pbs_file)
    with open(pbs_path, 'w+') as fh:
        fh.write(msg)
    return


def nan_check_misfits():
    msg = 'count=`grep "NaN" $RUNDIR/misfit*.dat | wc -l`\n'
    msg += 'if [ $count != 0 ]; then\n'
    msg += '        #stat=`grep -B 3 "NaN" $RUNDIR/misfit*.dat | '
    msg += "cut -d'-' -f1 | cut -d' ' -f1`\n"
    msg += '        comp=`grep -B1 "NaN" misfit*.dat '
    msg += "| cut -d'-' -f2 | "
    msg += 'awk "NR==1{print $1}"'
    msg += " | cut -d'.' -f1`\n"
    msg += '        eventid=`grep -B1 "NaN" misfit*.dat | '
    msg += "cut -d'-' -f3 | xargs`\n"
    msg += '        echo "misfit*.dat: NaN found in '
    msg += '$damp $iterno $comp $eventid" '
    msg += '>> $RUNDIR/progress.err\n'
    msg += '        echo "     NaN found, next damping" '
    msg += '>> $RUNDIR/progress.out\n'
    msg += '        exit 404\n'
    msg += 'fi\n\n'
    return msg


def dependencies_script(setup, uncertainties=False, iter_damp_test=False):
    """
    Creates a dependencies script to submit an inversion, using PBS arrays and
    the information given in the setup object.
    """

    if (
        setup.rundir.startswith('/quanta') or
        setup.rundir.startswith('/scratch') or
        setup.scripttype == 'slurm'
       ):
        msg = dorun_submission_header(setup, uncertainties=uncertainties)
        msg += 'cd $RUNDIR\n'
        if uncertainties is False:
            msg += 'rm inversion_out/inversion.db 2> /dev/null\n\n'

        msg += '#----------------#\n'
        msg += '#  DAMPING LOOP  #\n'
        msg += '#----------------#\n'
        if uncertainties is True:
            msg += 'damp=%s\n' % setup.damping
            msg += 'while read subset;\ndo\n'
        else:
            if iter_damp_test is False:
                msg += 'for damp in %s\n' % get_damping_string(setup)
            else:
                ds = get_damping_string(setup)
                ds = ds.split()[0] + ' #' + ' '.join(ds.split()[1:])
                msg += 'for damp in %s\n' % ds
            msg += 'do\n'
        if setup.intype == 'CST':
            msg += 'bash sbatch_startup.sh ${damp} $RUNDIR\n\n'
        elif setup.intype == 'FSI':
            msg += 'cd $RUNDIR\n'
            msg += dorun_preprundir(setup)
            msg += '\nN=0\n'
            sdirs = ''
            for segname, sdir in setup.segmentsdir.items():
                sdirs += "%s " % sdir

            msg += 'for sdir in %s; do\n' % sdirs

            msg += '    if (( $N >= 48 )); then\n'
            msg += '        for pid in ${pids[*]}; do\n'
            msg += '            wait $pid\n'
            msg += '            pids=${pids[@]/$pid}\n'
            msg += '                N=$(( N-1 ))\n'
            msg += '            break\n'
            msg += '        done\n'
            msg += '    fi\n'

            msg += '    bash sbatch_startup.sh $sdir &\n'

            msg += '    pids[${N}]=$!\n'
            msg += '    N=$(( N+1 ))\n'
            msg += 'done\n'
            msg += 'wait\n'

        msg += '#------------------------#\n'
        msg += '#  ITERATION LOOP START  #\n'
        msg += '#------------------------#\n'
        if iter_damp_test is False:
            msg += 'ed_it=%s\n' % (setup.iterations[1]+1)
        else:
            msg += 'ed_it=2 #%s\n' % (setup.iterations[1]+1)
        it = 'for (( iterno=%s; iterno<=${ed_it}; iterno++ ))\n'
        msg += it % (setup.iterations[0])
        msg += 'do\n\n'

        if len(setup.events) == 1:
            print('read_event_loop.pl will NOT fork correct!')
            print('read_event_loop.pl will mix up Set and cmt')
        msg += '#This part must be parallel\n'
        msg += 'echo "    -run events iteration $(( iterno-1 )) / '
        msg += '$(( ed_it-1 ))" >> $RUNDIR/progress.out\n'
        msg += 'perl read_event_loop.pl '
        msg += '$RUNDIR allevents 40 $iterno $ed_it $damp\n'
        msg += 'wait\n'
        # msg += 'if [ $? != 0 ]; then break; fi\n\n'

        msg += '\n# Wait until array finished, then do\n'
        if uncertainties is True:
            msg += 'bash uncertainties/sbatch_new.sh ${damp} ${iterno} '
            msg += '${ed_it} $subset $RUNDIR\n\n'
        else:
            msg += 'bash sbatch_new.sh ${damp} ${iterno} ${ed_it} $RUNDIR\n\n'

        if setup.rundir.startswith('/quanta'):
            msg += '\n# Sync after each iteration\n'
            if uncertainties is True:
                msg += 'cp -r $RUNDIR/uncertainties $MOM/.\n'
            else:
                msg += 'cp -r $RUNDIR/inversion_out $MOM/.\n'
            msg += 'cp $RUNDIR/*out $MOM/.\n'
            msg += 'cp $RUNDIR/*err $MOM/.\n'

        # msg += 'if [ $? != 0 ]; then break; fi\n\n'
        msg += '\n# Check here for NFS errors\n'
        msg += "err_no=`grep 'sqlite3.OperationalError'"
        if setup.rundir.startswith('/quanta'):
            msg += " $MOM/main*.err | wc -l`\n"
        else:
            msg += " $RUNDIR/main*.err | wc -l`\n"
        msg += "if [[ ! $err_no == 0 ]]; then\n"
        msg += "    echo 'NFS error, sqlite3 cannot write on disk' >> "
        msg += "$RUNDIR/progress.out\n"
        msg += "    exit 1\n"
        msg += "fi\n"
        msg += 'done\n'

        if uncertainties is True:
            msg += 'done < $RUNDIR/uncertainties/subsets\n'
        else:
            msg += 'done\n\n'
        msg += 'mv [0-9]*[A-Z].misfit inversion_out/.\n'
        msg += 'rm inversion.dat *misfit* mdcpl.out *.bin '
        msg += 'mnew* Cmat* Rmat* 2> /dev/null\n'
        if setup.keep_all_output is not True:
            msg += 'ls *dat |  xargs -I {{}} find {{}} -type f -not -name '
            msg += "'???????.dat' -delete \n"
        if setup.rundir.startswith('/scratch'):
            src = setup.rundir
            dst = setup.rundir.replace('/scratch', '/quanta1/home')
            dst = '/'.join(dst.split('/')[:-1])
            msg += 'SRC=%s\n' % (src)
            msg += 'DST=%s\n' % (dst)
            msg += 'mkdir -p %s\n' % dst
            msg += 'rsync -r --delete $SRC $DST'
        if setup.rundir.startswith('/quanta'):
            if uncertainties is True:
                msg += 'rsync -r --delete $RUNDIR/uncertainties $MOM/.\n'
            else:
                msg += 'rsync -r --delete $RUNDIR/inversion_out $MOM/.\n'
            msg += 'cp $RUNDIR/*out $MOM/.\n'
            msg += 'cp $RUNDIR/*err $MOM/.\n'
            msg += 'cd ~\n'
            msg += 'rm -rf /tmp/$USER\n'

    elif setup.intype == 'CST' and setup.rundir.startswith('/home'):
        start1 = "$(qsub -N %s -v damp=${damp} qsub_startup.pbs)\n"

        start2 = "$(qsub -N %s -v damp=${damp} -W depend=afterok:$OLD "
        start2 += "qsub_startup.pbs)\n"

        array = "$(qsub -N %s -v %s -W depend=afterok:$OLD "
        array += "qsub_array.pbs)\n"

        new = "$(qsub -N %s -v %s -W depend=afterokarray:$ARRAY "
        new += "qsub_new.pbs)\n"

        msg = '#!/bin/bash\n'
        msg += 'mkdir log\n'
        msg += 'st_it=%s\n' % setup.iterations[0]
        msg += 'ed_it=%s\n' % setup.iterations[1]
        msg += 'ed_it=$(($ed_it+1))\n'
        msg += 'count=1\n'

        var = 'damp=${damp},iterno=${iterno},ed_it=${ed_it}'

        msg += '\n'
        msg += '#----------------#\n'
        msg += '#  DAMPING LOOP  #\n'
        msg += '#----------------#\n'
        msg += 'for damp in %s\n' % get_damping_string(setup)
        msg += 'do\n'
        name = '%s_startup_d${damp}' % setup.pbsjob
        msg += 'if [[ count -eq 1 ]]; then\n'
        msg += 'START=' + start1 % name
        msg += 'else\n'
        msg += 'START=' + start2 % name

        msg += 'fi\n\n'
        msg += "STARTID=`echo $START | cut -d' ' -f4`\n"
        msg += 'OLD=$STARTID\n'
        msg += '((count++))\n'
        msg += '\n'
        msg += '#------------------------#\n'
        msg += '#  ITERATION LOOP START  #\n'
        msg += '#------------------------#\n'
        msg += 'for (( iterno=$st_it; iterno<=$ed_it; iterno++ ))\n'
        msg += 'do\n'

        name = '%s_array_d${damp}_it${iterno}' % setup.pbsjob
        msg += 'ARRAY=' + array % (name, var)
        msg += "ARRAYID=`echo $ARRAY| cut -d' ' -f4`\n"
        msg += 'ARRAY=$ARRAYID\n'
        name = '%s_new_d${damp}_it${iterno}' % setup.pbsjob
        msg += 'NEW=' + new % (name, var)
        msg += "NEWID=`echo $NEW | cut -d' ' -f4`\n"
        msg += 'OLD=$NEWID\n'
        msg += 'done\n'
        msg += 'done\n'

    elif setup.intype == 'FSI' and setup.rundir.startswith('/home'):
        msg = '#!/bin/bash\n'
        msg += 'mkdir log\n'
        msg += 'st_it=%s\n' % setup.iterations[0]
        msg += 'ed_it=%s\n' % setup.iterations[1]
        msg += 'ed_it=$(($ed_it+1))\n'
        msg += 'count=1\n'

        var = 'damp=${damp},iterno=${iterno},ed_it=${ed_it}'

        modes = ''
        for segname, sdir in setup.segmentsdir.items():
            if sdir.startswith('//nfs'):
                sdir = sdir.replace('//nfs/stig/', '/home/')
            modes = modes + "%s " % (sdir)

        msg += dorun_preprundir(setup)

        msg += '\n'
        msg += '#----------------#\n'
        msg += '#  DAMPING LOOP  #\n'
        msg += '#----------------#\n'
        msg += 'for damp in %s\n' % get_damping_string(setup)
        msg += 'do\n'
        msg += 'i=0\n'
        msg += 'for mdir in %s\n' % (modes)
        msg += 'do\n'
        msg += 'i=$(($i+1))\n'
        name = '%s_startup_d${damp}_i${i}' % setup.pbsjob
        msg += 'if [[ count -eq 1 ]]; then\n'
        msg += 'START=$(qsub -d $mdir -N %s -v damp=${damp} ' % name
        msg += 'qsub_startup.pbs)\n'
        msg += 'else\n'
        msg += 'START=$(qsub -d $mdir -N %s -v damp=${damp} ' % name
        msg += '-W depend=afterok:$OLD qsub_startup.pbs)\n'
        msg += 'fi\n'
        msg += 'done\n'
        msg += '((count++))\n'
        msg += 'OLD=$START\n'
        msg += 'while qstat -u jagt | grep -e "*startup_d*_i*" &> /dev/null\n'
        msg += 'do\n'
        msg += 'sleep 10\n'
        msg += 'done\n'
        msg += '\n'
        msg += '#------------------------#\n'
        msg += '#  ITERATION LOOP START  #\n'
        msg += '#------------------------#\n'
        msg += 'for (( iterno=$st_it; iterno<=$ed_it; iterno++ ))\n'
        msg += 'do\n'
        name = '%s_array_d${damp}_it${iterno}' % setup.pbsjob
        msg += 'ARRAY=$(qsub -N %s -v %s ' % (name, var)
        msg += '-W depend=afterok:$OLD qsub_array.pbs)\n'
        msg += '\n'

        name = '%s_new_d${damp}_it${iterno}' % setup.pbsjob
        msg += 'NEW=$(qsub -N %s -v %s ' % (name, var)
        msg += '-W depend=afterokarray:$ARRAY qsub_new.pbs)\n'
        msg += 'OLD=$NEW\n'
        msg += 'done\n'
        msg += 'done\n'

    if uncertainties is True:
        file_name = 'submit_uncertainties.sh'
    else:
        file_name = 'submit_dependencies.sh'

    file_path = join(setup.rundir, file_name)
    with open(file_path, 'w+') as fh:
        fh.write(msg)
    return


def pbsrundir_sanity_check(rundir):
    pbsrundir = rundir

    if pbsrundir.startswith('//nfs'):
        pbsrundir = pbsrundir.replace('//nfs/stig/', '/home/')
    #
    # if len(pbsrundir) > 45:
    #     pbsrundir = '"%s"\\\n"%s"' % (pbsrundir[0:46], pbsrundir[46:])
    return pbsrundir


def write_segment_files(mode, path=None, verbose=False):
    if path is None:
        path = '/quanta1/home/simons/splitting/modes/'
        path += mode

    comps = ['R', 'T', 'Z']

    for c in comps:
        if c == 'Z':
            seg = read_seg('db', channel=c, modes=format_name(mode).upper(),
                           min_snr=2)
        else:
            seg = read_seg('db', channel=c, modes=format_name(mode).upper(),
                           min_snr=1.5)

        if verbose:
            print(seg)
        os.chdir(os.path.join(path, "segments_%s" % c))
        seg.write('', eventfiles=True, overwrite=True, format='segment')

    return


def uncertainties_dir(setup, remove_existing=False, boot_strap=False,
                      jack_knive=True, size_of_subsets='full',
                      N_of_subsets=200, allevents_subset_ratio=None,
                      verbose=False):
    """
    Creates the uncertainties directory specifed for setup.rundir
    """
    if exists(join(setup.rundir, 'uncertainties')):
        if remove_existing is True:
            shutil.rmtree(join(setup.rundir, 'uncertainties'))
        else:
            raise IOError('uncertainties directory exists!')
    os.makedirs(join(setup.rundir, 'uncertainties'))

    src = join(setup.nmpydir, mergepath)
    dst = join(setup.rundir, "merge.py")
    shutil.copy(src, dst)

    with open(join(setup.rundir, 'allevents'), 'r') as fh:
        allevents = fh.readlines()

    file_subset_list = join(setup.rundir, 'uncertainties', 'subsets')
    subset_list = []

    N = len(allevents)
    if verbose:
        msg = '%s events read' % N
        print(msg)

    subsets = []
    wo_cmt = []
    cmt_list = ['031111B', '122604A', '022710A', '032805D', '062301E','060994A',
                '100494B', '081918A', '052413A', '031111Q', '122604Q','022710Q',
                '032805Q', '062301Q', '060994Q', '100494Q', '052413Q']
		        # '092503C', '081977B', '111506F', not needed

    if boot_strap is True:
        if size_of_subsets == 'full':
            size_of_subsets = N

        if size_of_subsets > N:
            size_of_subsets = N

        # Nrange = int(round(N/10.))
        N_cmt = 0

        for cmt in cmt_list:
            w_cmt = [x for x in allevents if cmt in x]
            if len(w_cmt) == 0:
                continue
            wo_cmt = set(allevents) - set(w_cmt)
            c = list(wo_cmt)
            c = [[item.rstrip() + ' ' + str(c.count(item))][0] for item in c]
            subsets += [c]
            N_cmt += 1

        for i in range(N_of_subsets-N_cmt):
            try:
                c = np.random.choice(wo_cmt, size=size_of_subsets,
                                     replace=False).tolist()
            except ValueError:
                if verbose is True:
                    print('size larger than allevents, replace set to True')
                c = np.random.choice(allevents, size=size_of_subsets,
                                     replace=True).tolist()
            c = [[item.rstrip() + ' ' + str(c.count(item))][0] for item in c]
            c = list(set(c))
            subsets += [c]

    elif jack_knive is True:
        # jack_knive should only remove one event at a time, nothing fancy
        allevents = [e.rstrip() for e in allevents]  # removving linebreak
        N_all = len(allevents)
        subsets = []

        if N_all < 20:
            for i, item in enumerate(allevents):
                subsets += [allevents[0:i] + allevents[i+1:]]
        else:
            # This loop removes creates subsets with one event only removed,
            # chosen from cmt_list
            for cmt in cmt_list:
                w_cmt = [x for x in allevents if cmt in x]
                wo_cmt = set(allevents) - set(w_cmt)
                subsets += [list(wo_cmt)]
            if verbose:
                msg = '%s fixed subsets created' % len(subsets)
                print(msg)

            if N_of_subsets > N_all:
                N = N_all
            else:
                N = N_of_subsets

            if allevents_subset_ratio is not None:
                size = int(allevents_subset_ratio * N_all)
            else:
                if size_of_subsets > N_all:
                    size = N_all - 1
                else:
                    size = size_of_subsets

            delta = int(abs(len(allevents) - size))
            if verbose:
                print('delta : %i' % delta)
            i = 0
            while(True):
                if i+delta > len(allevents):
                    break
                c = list(set(allevents) - set(allevents[i:i+delta]))
                subsets += [c]
                i= i + delta
            if verbose:
                msg = '%s subsets created in total' % len(subsets)
                print(msg)

    for i, subset in enumerate(subsets):
        fname = join(setup.rundir, 'uncertainties', "allevents_%s" % (i+1))
        # if verbose:
        #     print('subsets: %s' % fname)
        with open(fname, 'w') as fh:
            for line in subset:
                fh.write("%s\n" % line)
        with open(file_subset_list, 'a') as fh:
            subset_list += [fname]
            fh.write("%s\n" % fname)

    setup.uncertainties['files'] = subset_list

    return setup


def uncertainties_script(setup, damping):
    """
    creates:
        submission script:
        PATH/submit_uncertainties.sh

        new allevents files, containing the subsets:
        PATH/uncertainties/allevents_1-*

        List of the paths to these files:
        PATH/uncertainties/subsets

        modified sbatch_new.sh script:
        PATH/uncertainties/sbatch_new.sh

    """
    setup.damping = damping
    setup.uncertainties['damping'] = damping
    dependencies_script(setup, uncertainties=True)
    new_script(setup, uncertainties=True)
    array_script(setup, uncertainties=True)

    return setup
