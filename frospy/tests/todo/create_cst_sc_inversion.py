from frospy.scripts.inversion.create_input import create_std_sc_input
from frospy.core.setup.settings import Setup
from frospy.core.setup.builder import build_array_inversion
from frospy.scripts.inversion.submission import submit
import numpy as np

"""
param input['segmentsdir']: Segment-file or,
                            Path to folder that contains segments or,
                            Dictionary
"""

"""
Self coupling example
"""

"""
Sequential or array job
"""
mode = ("02T07", 6)
# modes = [("02T07", 14)]
modes_prem = None # ['00S19', '00T20', '01S13']  #, '10S04']  # , '00T24', '00S23']
# comps = [['R'], ['T'], ['Z']]

comp = ['R', 'T']
# startmodel = np.zeros(16)
# startmodel[0] = -8.
# startmodel[2] = 5.
# startmodel = 'S20RTS'
# startmodel = 'RR'
damp = np.logspace(-5, 1, 7).tolist()

# for mode, mp in zip(modes, modes_prem):
#     # for comp in comps:
#     if True:
#         print(mode, mp, comp)

# input = create_std_sc_input(mode, comp, startmodel, modes_prem=[mp],
#                             host='quanta')


for startmodel in ['PREM', 'S20RTS']:
    input = create_std_sc_input(mode, comp, startmodel, damp,
                                modes_prem=modes_prem,
                                host='quanta')
    if startmodel == 'S20RTS':
        input['rundir'] += '_S20'
        submit_run = False
        iter_damp_test = False
    else:
        submit_run = False
        iter_damp_test = False
        
    # input['rundir'] = '/quanta1/home/simons/splitting/modes/00t09/self_coupling/TR-Comp/18-0_plus'
    # input['bins'] = {'synseis_inv': 'synseis-invsyn',
    #                  'cst_partials': 'compu-deriv-cstCsyn'}
    setup = Setup('CST', input)
    # setup = build_inversion(input, 'CST', remove_existing=False, output='slurm')
    setup = build_array_inversion(input, 'CST', remove_existing=False,
                                  iter_damp_test=iter_damp_test)
    if submit_run is True:
        submit(setup.rundir)
