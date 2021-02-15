from frospy.core.spectrum.spectrum import Spectrum
import numpy as np
from frospy.util.base import update_progress
from frospy.util.read import read_st
from frospy.util.write import write_pickle

folder = '/tmp/eejit_simons/splitting/modes/01t08/self_coupling/'
folder += 'T-Comp/deg6/synthetics/'
data = folder + '010196C.ahx.syn'
syn = folder + 'broad/010196C.ahx.syn'

trace_ref = read_st(data)[0]
trace = read_st(syn)[0]


def timewindow_gridsearch(trace_ref, tw_ref, trace, fw):
    if len(trace_ref.data) != len(trace.data):
        print('file length does not match')
        return

    tend = trace_ref.stats.endtime - trace_ref.stats.starttime

    spec_ref = Spectrum(trace_ref, tw_ref[0], tw_ref[1])
    X_ref = spec_ref.data.fft.Data
    s_ref = spec_ref.flabel(fw[0])
    e_ref = spec_ref.flabel(fw[1])
    d2 = sum(abs(X_ref[s_ref:e_ref+1])**2.)

    eps = []
    # convert to hours
    loop_step = 3600.
    prog = 0
    progN = int(len(np.arange(0, tend, loop_step)) ** 2.) / 2
    progN += len(np.arange(0, tend, loop_step)) / 2
    for n in np.arange(0, tend, loop_step):
        prog += 1
        for m in np.arange(0, n, loop_step):
            prog += 1
            update_progress(float(prog)/progN)
            spec = Spectrum(trace, m/3600., n/3600.)
            X = spec.data.fft.Data
            s = spec.flabel(fw[0])
            e = spec.flabel(fw[1])
            cmf = sum(abs(X[s:e+1] - X_ref[s:e+1])**2.) / d2
            eps.append([m/3600., n/3600., cmf])

    return eps


eps = timewindow_gridsearch(trace_ref, [3, 60], trace, [2.270, 2.290])
write_pickle(eps, 'eps.pickle')


# For plotting

x = np.array(eps).reshape(len(eps), 3).transpose()
X = np.zeros((int(max(x[0]))+1, int(max(x[1]))+1))
for line in x.transpose():
    X[int(line[0])][int(line[1])] = line[2]
plt.imshow(X)
