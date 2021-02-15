from __future__ import absolute_import, print_function
from frospy.util.base import list2stream
import os

def sort_AttribDict(AD, key_type='int', reverse=False):
    """
    Sorts an AttribDict with Int or Str keys.
    """
    sorted_AD = []
    if key_type == 'int':
        keys = []
        for key in AD:
            keys.append(int(key))

        keys.sort(reverse=reverse)

        for key in keys:
            sorted_AD.append({str(key): AD[str(key)]})

    elif key_type == 'str':
        keys = []
        for key in AD:
            keys.append(key)
        keys.sort(reverse=reverse)

        for key in keys:
            sorted_AD.append({str(key): AD[key]})

    return sorted_AD


def sort_AttribDict_with_str_keys(AD, reverse=False):
    keys = []
    for key in AD:
        keys.append(int(key))

    keys.sort(reverse=reverse)
    sorted_AD = []

    for key in keys:
        sorted_AD.append({str(key): AD[str(key)]})

    return sorted_AD


def save_streamlist(streamlist, format='AH', filename=None, singlefiles=False):

    if singlefiles:
        for station in streamlist:
            time = station[0].stats.starttime
            name = station[0].stats.station
            network = station[0].stats.network
            try:
                location = station[0].stats.location
            except Exception:
                location = ''
            try:
                quality = station[0].stats.mseed['dataquality']
            except Exception:
                quality = ''

            fname = str(time.format_seed()).replace(",", ".") + "." + network \
                + "." + name + "." + location + "." + quality + "." + format

            station.write(fname, format=format)

    else:

        stream = list2stream(streamlist)
        stream.write(filename, format=format)


def write_cmt_file(cat, path=None):
    for event in cat:
        for ed in event.event_descriptions:
            if len(ed.text) == 7:
                cmt_id = ed.text

        origin = event.origins[0]
        lat = origin.latitude
        lon = origin.longitude
        depth_km = origin.depth / 1000.
        cmt = event.focal_mechanisms[0].moment_tensor.tensor

        m = '%.2f, %.2f, %.2f\n' % (lat, lon, depth_km)
        m += '%.3e, %.3e, %.3e, %.3e, %.3e, %.3e'.upper()
        m = m % (cmt.m_rr, cmt.m_tt, cmt.m_pp, cmt.m_rt, cmt.m_rp, cmt.m_tp)

        filename = "%s.cmt" % cmt_id
        if path is not None:
            filename = os.path.join(path, filename)

        with open(filename, 'w') as fh:
            fh.write(m)

    return
