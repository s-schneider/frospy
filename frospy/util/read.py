from __future__ import absolute_import
from future.utils import native_str
from obspy import read
# from frospy.core.modes import read as read_mode_class
from frospy.util.base import cat4stream, inv4stream
from frospy.util.ah import read_ahx
from frospy.util.array_util import (attach_network_to_traces,
                                  attach_coordinates_to_traces)
from obspy import read_events, read_inventory
from frospy import data as frospydata
import numpy as np
from obspy.core.util.attribdict import AttribDict
import obspy
import os
import struct
import pickle
import json


def read_st(input, format=None, inventory=None, event=None,
            net=None, client_name=None):

    try:
        stream = read(input, format)
    except UnicodeDecodeError:
        stream = read_ahx(input)

    try:
        if stream[0].stats._format == 'Q':
            for i, trace in enumerate(stream):
                stream[i].stats.distance = trace.stats.sh['DISTANCE']
                stream[i].stats.depth = trace.stats.sh['DEPTH']
                stream[i].stats.origin = trace.stats.sh['ORIGIN']
    except Exception:
        msg = 'No Q-file'

    if inventory:
        attach_network_to_traces(stream, inventory)
        if event:
            attach_coordinates_to_traces(stream, inventory, event)

    elif isinstance(net, str):
        if isinstance(client_name, str):
            try:
                cat = cat4stream(stream, client_name)
                inv = inv4stream(stream, net, client_name)

                attach_coordinates_to_traces(stream, inv, cat[0])
                attach_network_to_traces(stream, inv)
                print('Stream input with Meta-Information read.')
                return stream
            except Exception:
                msg = 'Error with Input: net or client_name wrong?'
                raise IOError(msg)

    return stream


def read_cmt(eventid=None):
    path = frospydata.__path__[0] + "/AD/cmt.dat"
    cmtlist = np.loadtxt(path, dtype={'names': ('id', 'cmt1', 'cmt2', 'cmt3',
                                                'cmt4', 'cmt5', 'cmt6',),
                                      'formats': ('|S15', np.float, np.float,
                                                  np.float, np.float, np.float,
                                                  np.float)},
                         delimiter=',').tolist()
    cmt = AttribDict()
    for i, entry in enumerate(cmtlist):
        tensor = AttribDict({
                            'Mrr': entry[1],
                            'Mtt': entry[2],
                            'Mpp': entry[3],
                            'Mrt': entry[4],
                            'Mrp': entry[5],
                            'Mtp': entry[6]
                            })
        cmt[entry[0]] = tensor
    if eventid:
        return cmt[eventid]
    else:
        return cmt


def read_cmt_file(file, full=False):
    with open(file) as fh:
        cmt = fh.readlines()
        for i, line in enumerate(cmt):
            cmt_tmp = line.strip().split()
            cmt[i] = [(float(j.split(',')[0])) for j in cmt_tmp]

    if full:
        return cmt
    else:
        return cmt[1]


def read_std_inv(inv='latlon'):
    """
    Reads stations stored in the file arwens_stations.xml . These stations
    were used for most of arwens work and represent a selection of high SNR
    stations.
    """
    if inv == 'full':
        inv_path = frospydata.__path__[0] + "/AD/arwens_stations.xml"
    else:
        inv_path = frospydata.__path__[0] + "/SAS/stations_latlon.xml"
    inv = read_inventory(inv_path)
    return inv


def read_std_cat(cmt_id=None):
    """
    Catalog file of 563 events.

    param cmt_id: string of event ID in cmt format
    """
    cat_path = frospydata.__path__[0] + "/AD/earthquakes_extended.xml"
    cat = obspy.Catalog()

    if type(cmt_id) == str and cmt_id.endswith('.cmt'):
        values = []
        with open(cmt_id, 'r') as fh:
            lines = fh.readlines()

        for i, l in enumerate(lines):
            # if csv values
            try:
                x = l.rstrip('\n').strip().split(',')
                values.append(np.array(x, dtype='float64'))
            # If tab separated values
            except ValueError:
                x = l.rstrip('\n').strip().split()
                values.append(np.array(x, dtype='float64'))

        event = obspy.core.event.Event()

        origin = obspy.core.event.Origin()
        origin.latitude = values[0][0]
        origin.longitude = values[0][1]
        origin.depth = values[0][2] * 1000.
        event.origins.append(origin)

        FM = obspy.core.event.FocalMechanism()
        mt = obspy.core.event.source.MomentTensor()
        mt.tensor = values[1][0]
        mt.tensor.m_tt = values[1][1]
        mt.tensor.m_pp = values[1][2]
        mt.tensor.m_rt = values[1][3]
        mt.tensor.m_rp = values[1][4]
        mt.tensor.m_tp = values[1][5]
        FM.moment_tensor = mt
        event.focal_mechanisms.append(FM)

        cat.append(event)

    else:
        cat_eq = read_events(cat_path)
        if cmt_id:
            for event in cat_eq:
                for desc in event.event_descriptions:
                    if desc.text in cmt_id:
                        cat.append(event)
                        break

            if len(cat) == 0:
                return None
        else:
            cat = cat_eq

    return cat


def read_dat_folder(folder):
    tmp_l = os.listdir(folder)
    flist = []

    for file in tmp_l:
        if file.endswith('.dat'):
            file = folder + file
            flist.append(file)
    flist.sort()

    return flist


def read_modes_in(modes_dir):
    """
    param modes_dir: path to directory containing:
                         modes.in
                         modes_cc.in
                     as defined for synseis hybrid

    type  modes_dir: string
    """
    if not modes_dir.endswith('/'):
        modes_dir += '/'

    fmodes_cc = modes_dir + 'modes_cc.in'
    fmodes = modes_dir + 'modes.in'

    try:
        with open(fmodes_cc, 'r') as fh:
            modes_cc = fh.readlines()

        if int(modes_cc[0]) == 0:
            modes_cc = [0]
        else:
            n = int(modes_cc[0])
            modes_cc = [x.strip() for x in modes_cc[0:n+1]]
    except IOError:
        modes_cc = None

    with open(fmodes, 'r') as fh:
        modes = fh.readlines()
    n = int(modes[0])  # number of modes in modes.in
    modes = [x.strip() for x in modes[0:n+1]]

    return modes, modes_cc


def write_mcst(cst, dst, modes_dir):
    # Read mode files
    modes, modes_cc = read_modes_in(modes_dir)
    sc_cdeg, sc_ddeg, cc_cdeg, cc_ddeg = get_mode_deg(modes, modes_cc)

    mcstfile = []
    for m, smax in zip(modes[1:], sc_cdeg):
        m = m.split()
        name = ''.join(m[:3]).upper()
        degs = cst[name].keys()

        for d in degs:
            if int(d) > int(smax):
                continue
            coeff = cst[name][d]
            for c in coeff:
                mcstfile.append(c)
        mcstfile.append(0)

    for m, smax in zip(modes_cc[1:], sc_cdeg):
        m = m.split()
        name1 = ''.join(m[:3]).upper()
        name2 = ''.join(m[3:6]).upper()
        name = '-'.join([name1, name2])
        degs = cst[name].keys()

        for d in degs:
            if int(d) > int(smax):
                continue
            coeff = cst[name][d]
            for c in coeff:
                mcstfile.append(c)
        mcstfile.append(0)

    N = len(mcstfile)
    with open('mcst_py.dat', 'w') as fh:
        fh.write("%s\n" % N)
        for line in mcstfile:
            fh.write("%s\n" % line)
    return


def read_omega_dat(qfile):
    omega = np.genfromtxt(qfile)
    return omega


def get_mode_names(modes, modes_cc):
    self_coupling = []
    cross_coupling = []
    for mode in modes[1:]:
        m = mode.split()
        self_coupling.append(''.join(m[0:3]))
    if modes_cc is None:
        cross_coupling = None
    else:
        for mode in modes_cc[1:]:
            m = mode.split()
            cross_coupling.append(''.join(m[0:3]) + '-' + ''.join(m[3:6]))
    return self_coupling, cross_coupling


def get_mode_deg(modes, modes_cc):
    # haydars format
    sc_cdeg = []
    sc_ddeg = []
    cc_cdeg = []
    cc_ddeg = []
    for mode in modes[1:]:
        m = mode.split()
        sc_cdeg.append(m[3])
        sc_ddeg.append(m[4])
    if modes_cc is None:
        cc_cdeg = cc_ddeg = None
    else:
        for mode in modes_cc[1:]:
            m = mode.split()
            cc_cdeg.append(m[-2])
            cc_ddeg.append(m[-1])
    return sc_cdeg, sc_ddeg, cc_cdeg, cc_ddeg


def get_mode_deg_dst(modes_sc_dst):
    # AD format without CC for dst
    sc_ddeg = []
    for mode in modes_sc_dst[1:]:
        m = mode.split()
        sc_ddeg.append(m[3])
    return sc_ddeg


def read_Total_VR(file):
    with open(file, 'r') as fh:
        content = fh.readlines()

    try:
        cmf = float(content[1].split()[-1])
        cmf = 1. - cmf/100.
    except ValueError:
        cmf = None
    try:
        amp_mf = float(content[2].split()[-1])
        amp_mf = 1. - amp_mf/100.
    except ValueError:
        amp_mf = None
    try:
        phase_mf = float(content[3].split()[-1])
        phase_mf = 1. - phase_mf/100
    except ValueError:
        phase_mf = None
    try:
        seg_cnt = int(content[4].split()[-1])
        npts = int(content[5].split()[-1])
    except ValueError:
        seg_cnt = None
        npts = None
    return cmf, amp_mf, phase_mf, seg_cnt, npts


def read_Total_MF(file):
    with open(file, 'r') as fh:
        content = fh.readlines()
    try:
        cmf = float(content[1].split()[-1])
    except ValueError:
        cmf = None
    try:
        amp_mf = float(content[2].split()[-1])
    except ValueError:
        amp_mf = None
    try:
        phase_mf = float(content[3].split()[-1])
    except ValueError:
        phase_mf = None
    try:
        seg_cnt = int(content[4].split()[-1])
        npts = int(content[5].split()[-1])
    except ValueError:
        seg_cnt = None
        npts = None
    return cmf, amp_mf, phase_mf, seg_cnt, npts


def read_varred_summary(file):
    with open(file, 'r') as fh:
        content = fh.readlines()

    cmf = float(content[1].split()[-1])
    cmf = 1. - cmf/100.
    return cmf


def read_atab_header_misfit(file, output='mf'):

    with open(file, mode='rb') as fhead:  # b is important -> binary
        fileContent = fhead.read()

    # mf = sum misfit^2
    misfit = struct.unpack('d', fileContent[4:12])[0]
    # dTd = sum data^2
    dTd = struct.unpack('d', fileContent[12:20])[0]

    if output == 'mf':
        av_mf = misfit/dTd
    elif output == 'vr':
        av_mf = (1. - misfit/dTd) * 100.

    return av_mf


def read_lines(file):
    with open(file, 'r') as fh:
        content = fh.readlines()
    for i, line in enumerate(content):
        content[i] = line.rstrip('\n')
    return content


def read_pickle(filename, **kwargs):
    """
    Read and return object from pickled object-file.

    :type filename: str
    :param filename: Name of the pickled object file to be read.
    :return: A python object.
    """
    kwargs = {}
    if isinstance(filename, (str, native_str)):  #, unicode)):
        with open(filename, 'rb') as fp:
            return pickle.load(fp, **kwargs)
    else:
        return pickle.load(filename, **kwargs)


def read_json(filename):
    if filename in ('AD', 'AD_cst', 'AD_cst.json'):
        path = "%s/AD/AD_cst.json" % frospydata.__path__[0]
    else:
        path = filename

    with open(path, 'r') as fh:
        data = json.load(fh)
    return data
