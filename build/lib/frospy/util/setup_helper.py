from collections import OrderedDict
from frospy.core.setup.settings import Setup


def create_setup_for_mode_smax(modes):
    '''
    type modes: list-like or OrderedDict
    param modes: list of lists, or tuples including the modes and
                 corresponding smax, e.g.:
        modes = (('0T4', 8), ('0S5', 10))
    '''
    if type(modes) in (list, tuple, OrderedDict):
        input_dict = {}
        if type(modes) == OrderedDict:
            input_dict['modes'] = modes
        elif type(modes) in (list, tuple):
            input_dict['modes'] = OrderedDict([modes])

        setup = Setup('CST', input_dict)
    return setup
