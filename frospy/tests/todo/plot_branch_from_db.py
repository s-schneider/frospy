from frospy.splitting.load import loadmodel
from frospy.postprocessing.plot import plot_coeffs_per_mode
from frospy.core.modes import format_name
from frospy.notebooks.data.SAS import get_paths
from frospy.core.splittingfunc import Set


dbpath = '/quanta1/home/simons/splitting/modes/cst.sqlite3'
paths = get_paths('all', results='all', host='quanta', dictionary=True)
SF = Set()


for mode, values in paths.items():
    # setup = Setup('CST', {'modes': {mode.upper(): values[1]}})
    damp = float(values[0].split('/')[-1].split('d')[1][:-1])
    try:
        SF += loadmodel(ifile=dbpath, modes=format_name(mode).upper(),
                        name='data', damp='0', db_model='data')
    except IndexError:
        continue

plot_coeffs_per_mode(SF_in=SF, cmap='grey', model=['RR', 'TZ'])

plot_coeffs_per_mode(SF_in=SF, fig_size=(12, 4), cmap='grey', degree=0,
                     plot_f_center=True, model=['RR', 'TZ', 'S20RTS'])
plot_coeffs_per_mode(SF_in=SF, fig_size=(12, 4), cmap='grey',
                     degree=0, plot_Q=True, model=['RR', 'TZ', 'AD'])

plot_coeffs_per_mode(SF_in=SF, fig_size=(18, 12), cmap='grey',
                     degree=2, model=['RR', 'TZ', 'S20RTS'])
