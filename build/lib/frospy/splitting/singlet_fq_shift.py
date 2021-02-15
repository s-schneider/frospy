#!/bin/python
"""
Module for visualizing shifts in frequency-Q of singlets,
resulting from using different models, coupling approximations, etc...

Required inputs:
- Two omega.dat files
- One modes.in file with the same modes as in omega.dat
"""
import numpy as np
import matplotlib.pyplot as plt
from frospy.core.modes import read as read_modes
from matplotlib.widgets import RectangleSelector
from matplotlib.gridspec import GridSpec
import os
import re
import subprocess


def plot_freq_Q(omega1=None, omega2=None, modefile=None, label1='omega1',
                label2='omega2', inspect=True):

    """
    example:
    plot_freq_Q('//nfs/stig/jagt/code/run_sc_syn0-3mHz.prem/Iter_1/matdiag/omega.dat',
                '//nfs/stig/jagt/code/run_syn0-3mHz.prem/Iter_1/matdiag/omega.dat',
                '//nfs/stig/jagt/code/run_syn0-3mHz.prem/config/modes.in',
                'Self coupling', 'Cross coupling')

    inspect provides a subplot with zoomed in version of main plot selection
    sens provides additionally a subplot with sensitivy kernel in selection
    """

    global freq, q, freq2, q2, pfreq, pq
    global ax, ax2, ax3, modename, fig
    global key3

    freq = []
    q = []
    freq2 = []
    q2 = []

    with open(omega1) as f:
        for line in f:
            cols = line.split()
            freq.append(float(cols[0]))
            q.append(1000./float(cols[1]))

    with open(omega2) as f:
        for line in f:
            cols = line.split()
            freq2.append(float(cols[0]))
            q2.append(1000./float(cols[1]))

    with open(modefile) as f:
        nmodes = int(f.readline().strip())
        modes = [next(f) for x in range(1, nmodes+1)]

    modelist = [mode.strip() for mode in modes]

    freq = np.array(freq)
    freq2 = np.array(freq2)
    q = np.array(q)
    q2 = np.array(q2)

    pfreq = []
    pq = []
    modename = []

    for mode in modelist:
        bits = mode.split(' ')
        lp = bits[2][-2:]
        if lp[0] == '0':
            L = lp[1:]
        else:
            L = lp
        nb = bits[0][-2:]
        if nb[0] == '0':
            N = nb[1:]
        else:
            N = nb

        name = str(N + bits[1] + L)

        try:
            mode_info = read_modes(mode_search=name)
            premfreq = mode_info[name.upper()]['freq']
            premq = mode_info[name.upper()]['Q']
            modename.append(name)

            pfreq.append(premfreq)
            pq.append(1000./premq)

        except KeyError:
            continue

    miny = min(q)
    maxy = max(q)
    minx = min(freq)
    maxx = max(freq)

    print(minx, maxx, miny, maxy)

    if inspect is True:
        gs = GridSpec(3, 3)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(gs[0:2, :])
        ax2 = fig.add_subplot(gs[2, 0:2])
        ax3 = fig.add_subplot(gs[2, 2])
    else:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    key1 = ax.scatter(freq, q, s=50, marker='o', facecolors='none',
                      edgecolors='cyan')
    key2 = ax.scatter(freq2, q2, s=10, marker='o', facecolors='none',
                      edgecolors='orange')
    key3 = ax.scatter(pfreq, pq, s=100, marker='o', facecolors='none',
                      edgecolors='magenta')

    for i in range(len(modename)):
        ax.text(pfreq[i], pq[i], s=modename[i], fontsize=8,
                verticalalignment='bottom')

    ax.set_title('Shift in singlet frequencies', fontsize=24)
    ax.set_xlabel('Frequency (mHz)')
    ax.set_ylabel('1000/Q')
    ax.legend([key1, key2, key3], [label1, label2, 'PREM (no coupling)'])

    if inspect is True:
        # set useblit True on gtkagg for enhanced performance
        span = RectangleSelector(ax, onselect, 'box', useblit=True,
                                 rectprops=dict(alpha=0.5, facecolor='red'))
    plt.show()


def onselect(eclick, erelease):

    # clear text and scatters in lower subplot
    # https://matplotlib.org/users/artists.html
    ax2.texts = []
    ax2.collections = []

    # selected box gives xmin,xmax,ymin,ymax
    xmin = eclick.xdata
    xmax = erelease.xdata
    ymin = eclick.ydata
    ymax = erelease.ydata

    # find modename of mode in selected box
    indspremf = np.where(np.logical_and(pfreq >= xmin, pfreq <= xmax))
    indspremq = np.where(np.logical_and(pq >= ymin, pq <= ymax))
    indprem = np.intersect1d(indspremf, indspremq)

    if len(indprem) == 1:
        moden = modename[int(indprem)]
        print(moden)
        ax2.legend([key3], [moden])
        plot_sens_kernel(moden)
        L = int(re.split('S|T', moden)[-1])
    elif len(indprem) < 1:
        print('No PREM modes found \n')
        ax2.legend([key3], ['No PREM modes in selected area'])
        ax3.lines = []
        L = 0
    elif len(indprem) > 1:
        print('Multiple prem modes found \n')
        ax2.legend([key3], ['Multiple modes in selected area'])
        ax3.lines = []
        L = 0

    # find indices of data in the selected box
    x = freq2
    y = q2
    indsx = np.where(np.logical_and(x >= xmin, x <= xmax))
    indsy = np.where(np.logical_and(y >= ymin, y <= ymax))
    inds = np.intersect1d(indsx, indsy)

    x2 = freq
    y2 = q
    indsx2 = np.where(np.logical_and(x2 >= xmin, x2 <= xmax))
    indsy2 = np.where(np.logical_and(y2 >= ymin, y2 <= ymax))
    inds2 = np.intersect1d(indsx2, indsy2)

    # check number of singlets against 2*l+1
    if len(inds) == 2*L+1 and len(indprem) == 1:
        print('Correct number of singlets \n')
    elif len(inds) != 2*L+1 and len(indprem) == 1:
        print('Incorrect number of singlets \n')

    # plotting
    thisx = x[inds]
    thisy = y[inds]
    thisx2 = x2[inds2]
    thisy2 = y2[inds2]
    ax2.scatter(thisx, thisy, s=10, marker='o', facecolors='none',
                edgecolors='orange')
    ax2.scatter(thisx2, thisy2, s=50, marker='o', facecolors='none',
                edgecolors='cyan')
    ax2.scatter(pfreq, pq, s=100, marker='o', facecolors='none',
                edgecolors='magenta')

    for i in indprem:
        ax2.text(pfreq[i], pq[i], s=modename[i], fontsize=10,
                 verticalalignment='bottom')

    try:
        ax2.set_xlim(min(thisx)-0.001, max(thisx)+0.001)
        ax2.set_ylim(min(thisy)-0.1, max(thisy)+0.1)
    except ValueError:
        print('No singlet data in selected area \n')

    fig.canvas.draw_idle()

    # save to textfile
    np.savetxt("text.out", np.c_[thisx, thisy])

    return


def plot_sens_kernel(moden):
    os.system('rm *-kernel.dat')

    ax3.lines = []

    bits = re.split('S|T', moden)

    if len(bits[0]) == 1:
        N = '00'+bits[0]
    elif len(bits[0]) == 2:
        N = '0'+bits[0]
    if len(bits[1]) == 1:
        L = '00'+bits[1]
    elif len(bits[1]) == 2:
        L = '0'+bits[1]

    if 'T' in moden:
        q = 'T'
    elif 'S' in moden:
        q = 'S'

    # get the kernel.dat files with modeplotqt
    os.system('echo "j" > input')
    os.system('echo %s %s %s >> input' % (N, q, L))
    res = subprocess.Popen('/net/home/jagt/modes/kernels/bin/modeplotqt <'
                           ' input', shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

    output, error = res.communicate()
    if res.returncode != 2:
        raise Exception(error)
        print('error: ', res.returncode)

    vp = []
    vs = []
    rho = []
    depth = []

    with open('rho-kernel.dat') as f:
        for line in f:
            cols = line.split()
            depth.append(float(cols[0]))
            rho.append(float(cols[1]))

    rhomax = np.amax(rho)

    with open('vs-kernel.dat') as f:
        for line in f:
            cols = line.split()
            vs.append(float(cols[1]))

    vsmax = np.amax(vs)

    intmax = np.amax([rhomax, vsmax])

    if q == 'S':
        with open('vp-kernel.dat') as f:
            for line in f:
                cols = line.split()
                vp.append(float(cols[1]))

        vpmax = np.amax(vp)
        intmax = np.amax([intmax, vpmax])

    ax3.set_ylim(0, 6400)
    ax3.set_xlim(-1.1*intmax, 1.1*intmax)

    plt.axhline(5700, color='k', lw=0.1)  # TZ
    plt.axhline(3480, color='k', lw=0.1)  # CMB
    plt.axhline(1220, color='k', lw=0.1)  # ICB
    plt.axvline(0, color='k', lw=0.1)

    han1, = ax3.plot(rho, depth, 'orange')
    han2, = ax3.plot(vs, depth, 'cyan')
    if q == 'S':
        han3, = ax3.plot(vp, depth, 'magenta')
        ax3.legend([han1, han2, han3], ['rho', 'vs', 'vp'])
    else:
        ax3.legend([han1, han2], ['rho', 'vs'])

    return
