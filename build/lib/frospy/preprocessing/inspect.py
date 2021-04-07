from __future__ import absolute_import, print_function
from frospy.util.base import (mask_data,
                              neighbouring_minima, get_local_extrema)


def remove_delta_pulses(trace, start_hour=10, delta_fac=5, verbose=False):
    """
    start_hour: start of timewindow to inspect trace in hours
    winsize_hours: timewindow size in hours
    """
    trace_work = trace.copy()
    start = int(start_hour * 3600. / trace_work.stats.delta)
    data = trace_work.data
    data_abs = abs(data[start:])

    a = get_local_extrema(data_abs, 'max')

    for x, y in zip(a, data_abs[a]):
        if y > data_abs.mean() * delta_fac:
            if verbose is True:
                print('removing outlier')
            mnm = neighbouring_minima(data_abs, x)
            s_ind = start + x + mnm[0]
            e_ind = start + x + mnm[1]
            data = mask_data(data, s_ind, e_ind + 1, 'linear')
    trace_work.data = data

    return trace_work
