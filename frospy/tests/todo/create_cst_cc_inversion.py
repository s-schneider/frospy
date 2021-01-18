from frospy.scripts.inversion.create_input import (create_std_cc_input,
                                                 get_modes_cc)
from frospy.core.setup.settings import Setup
from frospy.core.setup.builder import build_array_inversion
from obspy.core import AttribDict
import numpy as np


modes = [('00s10', 20), ('04S02', 4)]
startmodel = 'PREM'  # {'00S15': 'S20RTS', '01S11': 'S20RTS'}

# modes = [('00T07', 14), ('01T01', 2), ('00S07', 0), ('02S03', 0)]
# startmodel = 'PREM'

if len(modes) == 2:
    modes_cc = get_modes_cc(modes, [0, 0])
else:
    modes_cc = [
                ("%s-%s" % (modes[0][0], modes[1][0]), [0, 0]),
                ("%s-%s" % (modes[0][0], modes[2][0]), [0, 0]),
                ("%s-%s" % (modes[1][0], modes[2][0]), [0, 0])
                ]

damp = np.logspace(-5, 1, 7).tolist()
# damping = list(np.logspace(-5, 1, 7))[3:] + [0.025, 0.05, 0.075]
# damp = sorted(damping)

# --------------
# ----------PREM
# --------------
print('Building PREM run')
comps = ['Z']
# comps = ['R', 'T']
# comps = ['R', 'T', 'Z']

input = create_std_cc_input(modes, modes_cc, damp, comps=comps,
                            host='quanta', startmodel=startmodel)
                            # mode_components='all', host='quanta')
# job = 'c20c22'
job_pfix = '_fwAD'
input['rundir'] += job_pfix
# input['pbsjob'] += job_pfix
spath = '/quanta1/home/simons/splitting/modes/{}-{}/segments_Z_fwAD'
# spath = '/quanta1/home/simons/splitting/modes/{}-{}/segments_Z'
# spath = '/quanta1/home/simons/splitting/modes/{}-{}/segments_{}'

Zpath = spath.format(modes[0][0].lower(), modes[1][0].lower())
# spath = spath.format(modes[2][0].lower(), modes[0][0].lower())
# input['segmentsdir']['Z'] = spath


# Zpath = spath.format(modes[2][0].lower(), modes[0][0].lower(), 'Z')
# Rpath = spath.format(modes[2][0].lower(), modes[0][0].lower(), 'R')
# Tpath = spath.format(modes[2][0].lower(), modes[0][0].lower(), 'T')

# input['segmentsdir']['Z'] = Zpath
# input['segmentsdir']['R'] = Rpath
# input['segmentsdir']['T'] = Tpath

# spath = '/quanta1/home/simons/splitting/modes/{}-{}/segments_{}'
# for c in comps:
#     path = spath.format(modes[2][0].lower(), modes[0][0].lower(), c)
#     input['segmentsdir'][c] = path
input['segmentsdir']['Z'] = Zpath

setup = Setup('CST', input)
# setup.mzero['cst']['01T01'][0] = -4.
# setup.mzero['cst']['01T01'][-2] = 1.5
# setup.mzero['cst']['01T01'][-1] = -1.5
# setup.damping = [0.001]
# setup.iterations = [1, 10]
setup = build_array_inversion(setup, 'CST', remove_existing=False,
                              iter_damp_test=True)


# setup = Setup('CST', input)
# setup = build_array_inversion(setup, 'CST', remove_existing=False,
#                               iter_damp_test=False)
# mdamp = build_deg_dependent_damp(setup, modes=["00T16"], function=0.1, degs=[2, 4])
# damping_file(setup, mdamp)
