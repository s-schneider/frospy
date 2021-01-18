from __future__ import absolute_import, print_function
from frospy.core.modes import read as read_modes
from frospy.util.base import stream2array, taper_FT
from frospy.util.array_util import (attach_network_to_traces,
                                  attach_coordinates_to_traces,
                                  geometrical_center, find_closest_station,
                                  center_of_gravity)
from frospy.converter.bin2ascii import bin2ascii_matrix
from frospy.util.read import read_modes_in
from frospy.util.read import read_omega_dat


import numpy as np
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               ScalarFormatter)
from matplotlib.pyplot import cm
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
from itertools import cycle


from obspy.taup import TauPyModel
from obspy.core.event.event import Event
from obspy import Stream, Trace, Inventory
from obspy.core import AttribDict
from obspy.imaging.beachball import beach


class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here


def plot_freq(infile, outfile=None, newfigure=True, marker='D', size=10,
              fig_size=[7, 10], font_size=15,  fwin=False, label=None,
              limits=None, ticks=[None, None], title=None, papertype='A4',
              spacing=False, leg_loc=1, markmodes=True):
    """
    param	infile:	filename or list of filenames
    type	infile:	string or list of strings

    param limits: list of x and y limits
    type  limits: list-like, e.g. limits=[[1.855, 1.87], [321, 323]]

    example:

    ra25 = '/data/rad_ani/ra0.25/omega.dat'
    ra50 = '/data/rad_ani/ra0.5/omega.dat'
    ram25 = '/data/rad_ani/ra-0.25/omega.dat'
    ram50 = '/data/rad_ani/ra-0.5/omega.dat'
    sner = '/data/rad_ani/no_ra/omega.dat'
    ifile = [ra25, ra50, ram25, ram50, sner ]
    labels = ['rad-ani: 0.25', 'rad-ani: 0.5', 'rad-ani: -0.25',
              'rad-ani: -0.5', 'SNER \n+ rot \n+ ellip']
    size = [100,50,25,10,5]
    plot_freq(ifile, fwin=1., label=labels, limits=None, size=size,
              fig_size=False)
    """
    flen = len(infile)
    doplot = False
    if markmodes:
        modes = read_modes()
    else:
        modes = False

    if type(infile) == str:
        data = np.loadtxt(infile)
        doplot = True

    elif type(infile) == list and flen == 1:
        data = np.loadtxt(infile[0])
        doplot = True

    elif type(infile) == np.ndarray:
        doplot = True

    if doplot and not fwin:
        data.sort(axis=0)
        datat = infile.transpose()
        _doplot(datat, newfigure, marker, size, title, label, limits, leg_loc,
                markmodes, modes)

    elif doplot and fwin:
        data.sort(axis=0)
        datat = infile.transpose()
        start = 0.
        while True:
            end = start + fwin

            if end > data[0].max():
                end = data[0].max()

            fwindow = [[start, end], None]
            _doplot(datat, True, marker, size, title, label, fwindow, leg_loc,
                    markmodes, modes)
            start = end

            if start > data[0].max():
                break

    if type(infile) == list and flen != 1:
        if fwin:
            start = 0.
            # plt.ion()
            while True:
                fig, ax = plt.subplots()
                for i, file in enumerate(infile):
                    data = np.loadtxt(file)
                    datat = data.transpose()

                    end = start + fwin
                    if end > datat[0].max():
                        end = datat[0].max()
                    fwindow = [[start, end], None]

                    if type(marker) is not list:
                        marker_tmp = []
                        for _i in range(len(infile)):
                            marker_tmp.append(marker)
                        marker = marker_tmp

                    if type(size) is not list:
                        size_tmp = []
                        for _i in range(len(infile)):
                            size_tmp.append(size)
                        size = size_tmp

                    _doplot(datat, newfigure=False, marker=marker[i],
                            size=size[i], label=label[i], limits=fwindow,
                            markmodes=markmodes, modes=modes)
                start = end

                # plt.show()
                if start >= datat[0].max():
                    break
            # plt.ioff()
            # return

        else:
            for i, file in enumerate(infile):
                data = np.loadtxt(file)
                datat = data.transpose()
                if len(label) == len(infile):
                    if i == 0:
                        _doplot(datat, newfigure=False, marker=marker,
                                size=size, label=label[i], limits=limits,
                                markmodes=markmodes, modes=modes)
                    else:
                        _doplot(datat, newfigure=False, marker=marker,
                                size=size, label=label[i], limits=limits,
                                markmodes=False)
                else:
                    if i == 0:
                        _doplot(datat, newfigure=False, marker=marker,
                                size=size, label=label, limits=limits,
                                markmodes=markmodes)
                    else:
                        _doplot(datat, newfigure=False, marker=marker,
                                size=size, label=label, limits=limits,
                                markmodes=False)

    fig = plt.gcf()

    if spacing:
        fig.tight_layout(w_pad=.6, h_pad=.4)
    if fig_size:
        fig.set_size_inches(fig_size)

    mpl.rcParams.update({'font.size': font_size})
    if outfile:
        fig.savefig(outfile, dpi=400, orientation='portrait',
                    papertype=papertype)
        plt.close("all")

    else:
        plt.ion()
        print('No output format specified')
        plt.show()
        plt.ioff()


def _doplot(data, newfigure=True, marker='.', size=1, title='',
            label=None, limits=None, leg_loc=1, markmodes=False, modes=False):
        if newfigure:
            fig, ax = plt.subplots()
        else:
            ax = plt.gca()

        ax.set_title(title)
        ax.set_xlabel('Frequency (mHz)')
        ax.set_ylabel('Q')
        if limits:
            if limits[0]:
                ax.set_xlim(limits[0])
            if limits[1]:
                ax.set_ylim(limits[1])
        ax.scatter(data[0], data[1], s=size, label=label, marker=marker)
        ax.legend(loc=leg_loc)

        if markmodes and modes:
            trans = ax.get_xaxis_transform()
            for mode in modes:
                if mode.freq > limits[0][0] and mode.freq < limits[0][1]:
                        ax.annotate(r"$%s$" % mode.name,
                                    xy=(mode.freq, 1.04),
                                    xycoords=trans)
                        ax.axvline(x=mode.freq, linestyle='dashed',
                                   linewidth=1, color='grey')


def _domultiplot(data, x, y, gs, marker='.', size=1, title=None, label=None,
                 limits=None, ticks=None, newfigure=True):

    ax = plt.subplot(gs[x, y])
    ax.set_title(title)
    ax.set_ylabel('Q')
    ax.set_xlabel('Frequency (mHz)')

    if limits:
        if limits[0]:
            ax.set_xlim(limits[0])

        if ticks[0]:
            XmajorLocator = MultipleLocator(ticks[0][0])
            XmajorFormatter = FormatStrFormatter('%1.3f')
            XminorLocator = MultipleLocator(ticks[0][1])
            ax.xaxis.set_major_locator(XmajorLocator)
            ax.xaxis.set_major_formatter(XmajorFormatter)
            ax.xaxis.set_minor_locator(XminorLocator)

        if limits[1]:
            ax.set_ylim(limits[1])

        if ticks[1]:
            YmajorLocator = MultipleLocator(ticks[1][0])
            YmajorFormatter = FormatStrFormatter('%1.1f')
            YminorLocator = MultipleLocator(ticks[1][1])
            ax.yaxis.set_major_locator(YmajorLocator)
            ax.yaxis.set_major_formatter(YmajorFormatter)
            ax.yaxis.set_minor_locator(YminorLocator)

    plt.scatter(data[0], data[1], s=size, label=label, marker=marker)


def plot_single_spectrum(file, title=None, label=None, amplim=None,
                         linestyle=None, color=None, fontsize=None,
                         savefig=False, newfigure=True, xlim=None):
    """
    :param file: files to be plotted
    :type  file: list of strings
    :param label: list of labels to the corresponding files
    :type label: list of strings
    :param linestyle: list of linestyles to the corresponding files
    :type linestyle: list of strings
    :param color: list of colors to the corresponding files
    :type color: list of strings
    """
    if newfigure:
        fig, ax = plt.subplots()
    else:
        ax = plt.gca()
        fig = plt.gcf()

    if type(file) == str:
        file = [file]

    for i, fh in enumerate(file):
        data = np.loadtxt(fh).transpose()
        modes = read_modes()
        trans = ax.get_xaxis_transform()  # x in data untis, y in axes fraction

        if label and linestyle and color:
            ax.plot(data[0], data[1], label=label[i], color=color[i],
                    linestyle=linestyle[i][0], linewidth=linestyle[i][1])
        else:
            ax.plot(data[0], data[1])

    xmin, xmax = ax.get_xlim()
    for mode in modes:
        if mode.freq > xmin and mode.freq < xmax:
            if mode.name in [' _{0}S_{11} ']:
                ax.annotate(r"$%s$" % mode.name, xy=(mode.freq, 1.04),
                            xycoords=trans)
            else:
                ax.annotate(r"$%s$" % mode.name, xy=(mode.freq, 1.01),
                            xycoords=trans)
            ax.axvline(x=mode.freq, linestyle='dashed', linewidth=1,
                       color='grey')

    plt.legend(title=title)
    ax = format_exponent(ax)
    fig.set_size_inches(8, 7)
    ax.set_xlabel('Frequency (mHz)')
    ax.set_ylabel('Amplitude')
    mpl.rcParams.update({'font.size': fontsize})
    if xlim:
        ax.set_xlim(xlim)
    if amplim:
        ax.set_ylim(0, amplim)

    if savefig:
        fig.set_size_inches(10, 6)
        fig.savefig(savefig, dpi=300)
        plt.close("all")
    else:
        plt.ion()
        plt.show()
        plt.ioff()


def plot(st, inv=None, event=None, zoom=1, yinfo=False, stationlabel=True,
         epidistances=None, markphases=None, phaselabel=True,
         phaselabelclr='red', norm=False, clr=None, clrtrace=None,
         newfigure=True, savefig=False, dpi=400, xlabel=None, ylabel=None,
         yticks=False, t_axis=None, fs=15, tw=None, time_shift=None,
         verbose=False, kind='classic', labelfs=20, ylimit=None):
    """
    Alpha Version!

    Needs inventory and event of catalog for yinfo using

    param st: 	stream
    type  st:	obspy.core.stream.Stream

    param inv: inventory
    type  inv:

    param event: event of the seismogram
    type  event:

    param zoom: zoom factor of the traces
    type  zoom:	float

    param yinfo:	Plotting with y info as distance of traces
    type  yinfo:		bool

    param markphases: Phases, that should be marked in the plot, default is
                      "None"
    type  markphases: list

    param norm: Depiction of traces; unprocessed or normalized. Normalization
                options are:
                all - normalized on biggest value of all traces
                trace - each trace is normalized on its biggest value
    type  norm: string or bool

    param clr: Color of plot
    type  clr: string

    param clrtrace: dict containing tracenumber and color, e.g.:
                        >>>		clrtrace = {1: 'red', 7: 'green'}
                        >>>		trace 1 in red
                        >>>		trace 7 in green
                    or
                    attribute a trace, so only traces with that attribute are
                    plot in red, e.g.:
                        clrtrace='pocs', clrtrace='empty'

    type  clrtrace: list

    """

    # check for Data input
    if not isinstance(st, Stream):
        if not isinstance(st, Trace):
            try:
                if isinstance(yinfo, bool):
                    yinfo = 1
                plot_data(st, zoom=zoom, y_dist=yinfo, clr=clr,
                          newfigure=newfigure, savefig=savefig, dpi=dpi,
                          xlabel=xlabel, ylabel=ylabel, t_axis=t_axis, fs=fs)
                return
            except Exception:
                msg = "Wrong data input, must be Stream or Trace"
                raise TypeError(msg)
    if newfigure:
        # Set axis information and bools.
        fig, ax = plt.subplots()
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=fs)
        else:
            ax.set_xlabel("Time(s)", fontsize=fs)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=fs)

        ax.tick_params(axis='both', which='major', labelsize=fs)
        clr = 'black'
        cclr = 'black'

    else:
        ax = plt.gca()
        fig = plt.gcf()
        if not clr:
            clr = 'blue'
            cclr = 'blue'

    if isinstance(st, Stream):

        # Check if just specific timewindow should be plotted.
        try:
            tw = np.array(tw)
            twdelta = tw.max() - tw.min()
            t_axis = np.linspace(tw.min(), tw.max(),
                                 twdelta/float(st[0].stats.delta))
            npts_min = int(tw.min()/float(st[0].stats.delta))
            npts_max = int(tw.max()/float(st[0].stats.delta))
        except Exception:
            t_axis_max = st[0].stats.delta * st[0].stats.npts
            t_axis = np.linspace(0, t_axis_max, st[0].stats.npts)
            tw = np.array([0, t_axis_max])
            npts_min = 0
            npts_max = st[0].stats.npts

        data = stream2array(st)

        if time_shift:
            data = np.roll(data, time_shift*int(st[0].stats.sampling_rate))

        spacing = 2.
        ax.set_xlim(tw.min(), tw.max())
        isinv = False
        isevent = False

        # Check if inventory and catalog is input,
        # then calculate distances etc.
        if isinstance(inv, Inventory):
            isinv = True
        if isinstance(event, Event):
            isevent = True

        if isinv:
            # Calculates y-axis info using epidistance information of stream.
            # Check if there is a network entry
            attach_network_to_traces(st, inv)
            attach_coordinates_to_traces(st, inv, event)

        else:
            for trace in st:
                try:
                    if not trace.stats.distance:
                        isinv = False
                        break
                    else:
                        isinv = True
                except Exception:
                        isinv = False
                        break

        try:
            depth = event.origins[0]['depth']/1000.
            isevent = True
        except AttributeError:
            try:
                depth = st[0].stats.depth
                isevent = True
            except AttributeError:
                isevent = False

        yold = 0

        # Normalize Data, if set to 'all'
        if norm in ['all']:
            data = data/data.max()

        if yinfo:
            ymin = st[0].stats.distance
            ymax = st[0].stats.distance

        if markphases and isinv and isevent:
            try:
                origin = event.origins[0]['time']
                depth = event.origins[0]['depth']/1000.
            except Exception:
                origin = st[0].stats.origin
                depth = st[0].stats.depth

            m = TauPyModel('ak135')

        if kind in ('classic', 'Classic'):
            N = len(data)
            for j, trace in enumerate(data):

                # Normalize trace, if set to 'trace'
                if norm in ['trace']:
                    trace = trace/trace.max()

                try:
                    y_dist = st[j].stats.distance
                except Exception:
                    y_dist = yold + 1

                if markphases and isinv and isevent:

                    arrivals = m.get_travel_times(depth, y_dist,
                                                  phase_list=markphases)
                    timetable = [[], []]
                    for k, phase in enumerate(arrivals):
                        phase_name = phase.name
                        t = phase.time
                        phase_time = origin + t - st[j].stats.starttime
                        Phase_npt = int(phase_time/st[j].stats.delta)
                        Phase = Phase_npt * st[j].stats.delta

                        if Phase < t_axis.min() or Phase > t_axis.max():
                            continue
                        else:
                            timetable[0].append(phase_name)
                            timetable[1].append(Phase)

                    if not timetable[0] or not timetable[1]:
                        print('Phases not in Seismogram')
                        plt.close('all')
                        return

                    if yinfo:
                        if not ylabel:
                            ax.set_ylabel("Distance(deg)", fontsize=fs)

                        if st[j].stats.distance < ymin:
                            ymin = st[j].stats.distance
                        if st[j].stats.distance > ymax:
                            ymax = st[j].stats.distance

                        try:
                            if j in clrtrace:
                                cclr = clrtrace[j]
                            else:
                                cclr = clr
                        except TypeError:
                            if clrtrace in st[j].stats:
                                cclr = 'red'
                            else:
                                cclr = clr
                        except Exception:
                            cclr = clr

                        if stationlabel:
                            ax.annotate('%s' % st[j].stats.station,
                                        xy=(1 + tw.min(), y_dist+0.1))

                        ax.plot(t_axis, zoom*trace[npts_min: npts_max]+y_dist,
                                color=cclr)
                        ax.plot((timetable[1], timetable[1]),
                                (-1+y_dist, 1+y_dist), color=phaselabelclr)
                        if verbose:
                            print(timetable[1] + st[j].stats.shifttime)
                            ax.plot((timetable[1] + st[j].stats.shifttime,
                                    timetable[1] + st[j].stats.shifttime),
                                    (-1+y_dist, 1+y_dist), color=phaselabelclr)

                        if phaselabel:
                            for time, key in enumerate(timetable[0]):
                                ax.annotate('%s' % key,
                                            xy=(timetable[1][time], y_dist))
                                if verbose:
                                    tt = timetable[1][time]
                                    tt += st[j].stats.shifttime
                                    ax.annotate('%s' % key, xy=(tt, y_dist))
                        else:
                            continue

                    else:
                        if not ylabel:
                            ax.set_ylabel("No. of trace", fontsize=fs)

                        try:
                            if j in clrtrace:
                                cclr = clrtrace[j]
                            else:
                                cclr = clr
                        except TypeError:
                            if clrtrace in st[j].stats:
                                cclr = 'red'
                            else:
                                cclr = clr
                        except Exception:
                            cclr = clr

                        fig.gca().yaxis.set_major_locator(plt.NullLocator())

                        if stationlabel:
                            ax.annotate('%s' % st[j].stats.station,
                                        xy=(1 + tw.min(), spacing*j+0.1))

                        ax.plot(t_axis,
                                zoom*trace[npts_min: npts_max]+spacing*j,
                                color=cclr)
                        ax.plot((timetable[1], timetable[1]),
                                (-1+spacing*j, 1+spacing*j),
                                color=phaselabelclr)

                        if verbose:
                            print(st[j].stats.shifttime)
                            print(timetable[1] + st[j].stats.shifttime)

                            ax.plot((timetable[1] + st[j].stats.shifttime,
                                    timetable[1] + st[j].stats.shifttime),
                                    (-1+spacing*j, 1+spacing*j),
                                    color=phaselabelclr)

                        if phaselabel:
                            for time, key in enumerate(timetable[0]):
                                ax.annotate('%s' % key,
                                            xy=(timetable[1][time], spacing*j))
                                if verbose:
                                    tt = timetable[1][time]
                                    tt += st[j].stats.shifttime
                                    ax.annotate('%s' % key,
                                                xy=(tt, spacing*j))
                        else:
                            continue

                elif markphases and not isinv:
                    msg = 'markphases needs Inventory Information, not found.'
                    raise IOError(msg)

                elif markphases and not isevent:
                    msg = 'markphases needs Event Information, not found.'
                    raise IOError(msg)

                elif (type(epidistances) == np.ndarray or
                      type(epidistances) == list):
                    y_dist = epidistances
                    if not ylabel:
                        ax.set_ylabel("Distance(deg)", fontsize=fs)
                    try:
                        if j in clrtrace:
                            cclr = clrtrace[j]
                        else:
                            cclr = clr
                    except Exception:
                        cclr = clr

                    if stationlabel:
                        ax.annotate('%s' % st[j].stats.station,
                                    xy=(1 + tw.min(), y_dist[j]+0.1))

                    ax.plot(t_axis, zoom*trace[npts_min: npts_max] + y_dist[j],
                            color=cclr)

                else:
                    if yinfo:

                        try:
                            if not ylabel:
                                ax.set_ylabel("Distance(deg)", fontsize=fs)

                            try:
                                if j in clrtrace:
                                    cclr = clrtrace[j]
                                else:
                                    cclr = clr
                            except TypeError:
                                if clrtrace in st[j].stats:
                                    cclr = 'red'
                                else:
                                    cclr = clr
                            except Exception:
                                cclr = clr

                        except Exception:
                            msg = 'Oops, something not found.'
                            raise IOError(msg)

                        if st[j].stats.distance < ymin:
                            ymin = st[j].stats.distance
                        if st[j].stats.distance > ymax:
                            ymax = st[j].stats.distance

                        if stationlabel:
                            ax.annotate('%s' % st[j].stats.station,
                                        xy=(1 + tw.min(), y_dist+0.1))

                        ax.plot(t_axis, zoom * trace[npts_min: npts_max] +
                                y_dist, color=cclr)

                    else:
                        if not ylabel:
                            ax.set_ylabel("No. of trace", fontsize=fs)
                        try:
                            if j in clrtrace:
                                cclr = clrtrace[j]
                            else:
                                cclr = clr
                        except Exception:
                            cclr = clr

                        if stationlabel:
                            xran = abs(ax.get_xlim()[0] - ax.get_xlim()[1])
                            trans = ax.get_xaxis_transform()
                            ax.annotate('%s' % st[j].stats.station,
                                        xy=(-0.08*xran, j/N),
                                        xycoords=trans)

                        ax.plot(t_axis, zoom * trace[npts_min: npts_max] +
                                spacing * j, color=cclr)

                yold = y_dist

            if not yticks:
                fig.gca().yaxis.set_major_locator(plt.NullLocator())

        elif kind in ['contour', 'Contour']:
            yrange = np.arange(len(st))

            try:
                cax = ax.imshow(zoom*data[tw.min():tw.max()], aspect='auto',
                                extent=(tw.min(), tw.max(), yrange.min(),
                                yrange.max()), origin='lower', cmap='seismic')
            except Exception:
                cax = ax.imshow(zoom*data, aspect='auto', extent=(t_axis.min(),
                                t_axis.max(), yrange.min(), yrange.max()),
                                origin='lower', cmap='seismic')

            if markphases and isinv and isevent:

                stats = np.arange(len(st))
                for j in stats:

                    arrivals = m.get_travel_times(depth, st[j].stats.distance,
                                                  phase_list=markphases)
                    timetable = [[], []]

                    for phase in arrivals:
                        t = phase.time
                        phase_time = origin + t - st[j].stats.starttime
                        Phase_npt = int(phase_time / st[j].stats.delta)
                        tPhase = Phase_npt * st[j].stats.delta
                        name = phase.name

                        if tPhase > t_axis.max() or tPhase < t_axis.min():
                            continue

                        ax.autoscale(False)

                        ax.plot((tPhase, tPhase),
                                (-labelfs/100.+j, labelfs/100.+j),
                                color='red')

                        if j == stats.max():
                            if name in [u'P^220P']:
                                ax.annotate('$P^{220}P$', xy=(tPhase, j),
                                            xytext=(tPhase + 1, j+1.2),
                                            fontsize=labelfs, color='black')

                            elif name in [u'P^410P']:
                                ax.annotate('$P^{410}P$', xy=(tPhase, j),
                                            xytext=(tPhase + 1, j+1.2),
                                            fontsize=labelfs, color='black')

                            elif name in [u'P^660P']:
                                ax.annotate('$P^{660}P$', xy=(tPhase, j),
                                            xytext=(tPhase + 1, j+1.2),
                                            fontsize=labelfs, color='black')

                            else:
                                ax.annotate('$%s$' % name, xy=(tPhase, j),
                                            xytext=(tPhase + 1, j+1.2),
                                            fontsize=labelfs, color='black')

            cbar = fig.colorbar(cax, format='%.1f')
            cbar.set_clim(-1, 1)
            cbar.ax.tick_params(labelsize=fs)
            cbar.ax.set_ylabel('A', fontsize=fs)

        if yinfo:
            ylim = (ymin-1, ymax+1)
            ax.set_ylim(ylim)

        if savefig:
            fig.set_size_inches(8, 7)
            fig.savefig(savefig, dpi=dpi)
            plt.close("all")
        else:
            plt.ion()
            plt.draw()
            plt.show()
            plt.ioff()

    elif isinstance(st, Trace):

        t_axis = np.linspace(0, st.stats.delta * st.stats.npts, st.stats.npts)
        data = st.data.copy()

        if norm in ['all', 'All', 'trace', 'Trace']:
            data = data/data.max()
        try:
            y_dist = st.stats.distance
        except Exception:
            msg = "No distance information attached to trace"
            msg += ", no phases are calculated!"
            print(msg)
            markphases = False

        if ylimit:
            ax.set_ylim(ylimit)

        if markphases:
            try:
                origin = event.origins[0]['time']
                depth = event.origins[0]['depth']/1000.
            except Exception:
                origin = st.stats.origin
                depth = st.stats.depth

            m = TauPyModel('ak135')
            arrivals = m.get_travel_times(depth, y_dist, phase_list=markphases)
            timetable = [[], []]
            for k, phase in enumerate(arrivals):
                phase_name = phase.name
                t = phase.time
                phase_time = origin + t - st.stats.starttime
                Phase_npt = int(phase_time/st.stats.delta)
                Phase = Phase_npt * st.stats.delta

                if Phase < t_axis.min() or Phase > t_axis.max():
                    continue
                else:
                    timetable[0].append(phase_name)
                    timetable[1].append(Phase)

            if not ylabel:
                ax.set_ylabel("Amplitude", fontsize=fs)
            title = st.stats.network+'.'+st.stats.station+'.'
            title += st.stats.location+'.'+st.stats.channel

            ax.set_title(title, fontsize=fs)
            ax.plot(t_axis, zoom*data, color=clr)
            ax.plot((timetable[1], timetable[1]), (-0.5, 0.5),
                    color=phaselabelclr)
            if phaselabel:
                for time, key in enumerate(timetable[0]):
                    ax.annotate('%s' % key, xy=(timetable[1][time]+5, 0.55))

        else:
            if not ylabel:
                ax.set_ylabel("Amplitude", fontsize=fs)
            title = st.stats.network+'.'+st.stats.station+'.'
            title += st.stats.location+'.'+st.stats.channel
            ax.set_title(title, fontsize=fs)
            ax.plot(t_axis, zoom*data, color=clr)

        if savefig:
            plt.savefig(savefig)
            plt.close("all")
        else:
            plt.ion()
            plt.draw()
            plt.show()
            plt.ioff()


def plot_data(data, zoom=1, y_dist=1, label=None, clr=None, newfigure=True,
              savefig=False, dpi=400, xlabel=None, ylabel=None, t_axis=None,
              fs=15):
    """
    Alpha Version!
    Time axis has no time-ticks --> Working on right now

    Needs inventory and catalog for yinfo using
    param st: 	array of data
    type st:	np.array
    param zoom: zoom factor of the traces
    type zoom:	float
    param y_dist:	separating distance between traces, for example equidistant
                    with "1"
                    or import epidist-list via epidist
    type y_dist:	int or list
    """
    if newfigure:
        fig, ax = plt.subplots()
        ax.set_xlabel(xlabel, fontsize=fs)
        ax.set_ylabel(ylabel, fontsize=fs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
        ticks = mpl.ticker.FuncFormatter(lambda r, pos: '{0:g}'.format(r /
                                                                       y_dist))
        ax.yaxis.set_major_formatter(ticks)
        clr = 'black'

    else:
        ax = plt.gca()
        fig = plt.gcf()
        if not clr:
            clr = 'blue'

    if data.ndim == 1:
        if not t_axis:
            t_axis = np.arange(len(data))
        ax.plot(t_axis, zoom*data, color=clr, label=label)
    else:
        if not t_axis:
            t_axis = np.arange(len(data[0]))

        for i, trace in enumerate(data):
            if isinstance(y_dist, int):
                try:
                    if i == 0:
                        ax.plot(t_axis, zoom*trace + y_dist*i, color=clr,
                                label=label)
                    else:
                        ax.plot(t_axis, zoom*trace + y_dist*i, color=clr)
                except Exception:
                    if i == 0:
                        ax.plot(zoom*trace + y_dist*i, color=clr, label=label)
                    else:
                        ax.plot(zoom*trace + y_dist*i, color=clr)

    if savefig:
        fig.set_size_inches(12, 10)
        fig.savefig(savefig, dpi=dpi)
        plt.close('all')
    else:
        plt.ion()
        plt.draw()
        ax.legend()
        plt.show()
        plt.ioff()


def beachball(fm, linewidth=2, facecolor='b', bgcolor='w', edgecolor='k',
              alpha=1.0, xy=(0, 0), width=200, size=100, nofill=False,
              zorder=100, axis=None):
        """
        Draws a beach ball diagram of an earthquake focal mechanism.

        S1, D1, and R1, the strike, dip and rake of one of the focal planes,
        can be vectors of multiple focal mechanisms.

        :param fm: Focal mechanism that is either number of mechanisms (NM) by
        3 (strike, dip, and rake) or NM x 6 (M11, M22, M33, M12, M13, M23 - the
        six independent components of the moment tensor, where the coordinate
        system is 1,2,3 = Up,South,East which equals r,theta,phi). The strike
        is of the first plane, clockwise relative to north.
        The dip is of the first plane, defined clockwise and perpendicular to
        strike, relative to horizontal such that 0 is horizontal and 90 is
        vertical. The rake is of the first focal plane solution. 90 moves the
        hanging wall up-dip (thrust), 0 moves it in the strike direction
        (left-lateral), -90 moves it down-dip (normal), and 180 moves it
        opposite to strike (right-lateral).
        :param facecolor: Color to use for quadrants of tension; can be a
            string, e.g. ``'r'``, ``'b'``
            or three component color vector, [R G B].
            Defaults to ``'b'`` (blue).
        :param bgcolor: The background color. Defaults to ``'w'`` (white).
        :param edgecolor: Color of the edges. Defaults to ``'k'`` (black).
        :param alpha: The alpha level of the beach ball. Defaults to ``1.0``
            (opaque).
        :param xy: Origin position of the beach ball as tuple. Defaults to
            ``(0, 0)``.
        :type width: int
        :param width: Symbol size of beach ball. Defaults to ``200``.
        :param size: Controls the number of interpolation points for the
            curves. Minimum is automatically set to ``100``.
        :param nofill: Do not fill the beach ball, but only plot the planes.
        :param zorder: Set zorder. Artists with lower zorder values are drawn
            first.
        :param outfile: Output file string. Also used to automatically
            determine the output format. Supported file formats depend on your
            matplotlib backend. Most backends support png, pdf, ps, eps and
            svg. Defaults to ``None``.
        :param format: Format of the graph picture. If no format is given the
            outfile parameter will be used to try to automatically determine
            the output format. If no format is found it defaults to png output.
            If no outfile is specified but a format is, than a binary
            imagestring will be returned.
            Defaults to ``None``.
        :param ax: Give an existing ax instance to plot into. New Figure if
            set to ``None``.
        """
        plot_width = width * 0.95

        # plot the figure
        if not axis:
            fig = plt.figure(figsize=(3, 3), dpi=100)
            fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
            fig.set_figheight(width // 100)
            fig.set_figwidth(width // 100)
            ax = fig.add_subplot(111, aspect='equal')
        else:
            ax = axis

        # hide axes + ticks
        ax.axison = False

        # plot the collection
        collection = beach(fm, linewidth=linewidth, facecolor=facecolor,
                           edgecolor=edgecolor, bgcolor=bgcolor,
                           alpha=alpha, nofill=nofill, xy=xy,
                           width=plot_width, size=size, zorder=zorder)
        ax.add_collection(collection)

        ax.autoscale_view(tight=True, scalex=True, scaley=True)
        ax.patch.set_facecolor('red')
        ax.patch.set_alpha(0.0)

        if axis:
            return
        else:
            return fig, ax


def axvlines(xs, ax, **plot_kwargs):
    """
    Use this for modes plot
    Draw vertical lines on plot
    :param xs: A scalar, list, or 1D array of horizontal offsets
    :param plot_kwargs: Keyword arguments to be passed to plot
    :return: The plot object corresponding to the lines.
    """
    xs = np.array((xs, ) if np.isscalar(xs) else xs, copy=False)
    lims = ax.get_ylim()
    x_points = np.repeat(xs[:, None], repeats=3, axis=1).flatten()
    y_points = np.repeat(np.array(lims + (np.nan, ))[None, :], repeats=len(xs),
                         axis=0).flatten()

    if ax:
        plot = ax.plot(x_points, y_points, scaley=False, **plot_kwargs)
    else:
        plot = plt.plot(x_points, y_points, scaley=False, **plot_kwargs)
    return plot


def axtexts(xs, texts, ax, **plot_kwargs):
    """
    Use this for modes plot
    Draw vertical lines on plot
    :param xs: A scalar, list, or 1D array of horizontal offsets
    :param plot_kwargs: Keyword arguments to be passed to plot
    :return: The plot object corresponding to the lines.
    """

    for x, text in zip(xs, texts):
        if ax:
            plot = ax.text(x, s=text, **plot_kwargs)
        else:
            plot = plt.text(x, s=text, **plot_kwargs)

    return plot


def axhlines(ys, ax, **plot_kwargs):
    """
    Draw horizontal lines across plot
    :param ys: A scalar, list, or 1D array of vertical offsets
    :param plot_kwargs: Keyword arguments to be passed to plot
    :return: The plot object corresponding to the lines.
    """
    ys = np.array((ys, ) if np.isscalar(ys) else ys, copy=False)
    lims = plt.gca().get_xlim()
    y_points = np.repeat(ys[:, None], repeats=3, axis=1).flatten()
    x_points = np.repeat(np.array(lims + (np.nan, ))[None, :], repeats=len(ys),
                         axis=0).flatten()
    if ax:
        plot = ax.plot(x_points, y_points, scalex=False, **plot_kwargs)
    else:
        plot = plt.plot(x_points, y_points, scalex=False, **plot_kwargs)
    return plot


def format_exponent(ax, axis='y'):

    # Change the ticklabel format to scientific format
    # ax.ticklabel_format(axis=axis, style='sci', scilimits=(-2, 2))
    ScalarFormatterForceFormat()
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((-2,2))
    ax.yaxis.set_major_formatter(yfmt)

    # Get the appropriate axis
    if axis == 'y':
        ax_axis = ax.yaxis
        x_pos = 0.0
        y_pos = 1.0
        horizontalalignment = 'left'
        verticalalignment = 'bottom'
    else:
        ax_axis = ax.xaxis
        x_pos = 1.0
        y_pos = -0.05
        horizontalalignment = 'right'
        verticalalignment = 'top'

    # Run plt.tight_layout() because otherwise the offset text doesn't update
    plt.tight_layout()
    # THIS IS A BUG
    # Well, at least it's sub-optimal because you might not
    # want to use tight_layout(). If anyone has a better way of
    # ensuring the offset text is updated appropriately
    # please comment!

    # Get the offset value
    offset = ax_axis.get_offset_text().get_text()

    if len(offset) > 0:
        # Get that exponent value and change it into latex format
        minus_sign = u'\u2212'
        expo = np.float(offset.replace(minus_sign, '-').split('e')[-1])
        offset_text = r'x$\mathregular{10^{%d}}$' % expo

        # Turn off the offset text that's calculated automatically
        ax_axis.offsetText.set_visible(False)
        ax.text(x_pos, y_pos, offset_text, transform=ax.transAxes,
                horizontalalignment=horizontalalignment,
                verticalalignment=verticalalignment)
    return ax


def plot_inv(inventory, projection="local"):
    """
    Function to plot the geometry of the array,
    including its center of gravity and geometrical center

    :type inventory: obspy.core.inventory.inventory.Inventory
    :param inventory: Inventory to be plotted

    :type projection: strg, optional
    :param projection: The map projection. Currently supported are:

    * ``"global"`` (Will plot the whole world.)
    * ``"ortho"`` (Will center around the mean lat/long.)
    * ``"local"`` (Will plot around local events)
    """
    if isinstance(inventory, Inventory):
        inventory.plot(projection, show=False)
        bmap = plt.gca().basemap

        grav = center_of_gravity(inventory)
        x, y = bmap(grav["longitude"], grav["latitude"])
        bmap.scatter(x, y, marker="x", c="red", s=40, zorder=20)
        plt.text(x, y, "Center of Gravity", color="red")

        geo = geometrical_center(inventory)
        x, y = bmap(geo["longitude"], geo["latitude"])
        bmap.scatter(x, y, marker="x", c="green", s=40, zorder=20)
        plt.text(x, y, "Geometrical Center", color="green")
        plt.ion()
        plt.draw()
        plt.show()
        plt.ioff()


def plot_vespa(data, st, inv=None, event=None,
               markphases=['ttall', 'P^410P', 'P^660P'], plot='classic',
               cmap='seismic', tw=None, savefig=False, dpi=400, fs=25, power=4,
               marker='|', markerclr='red', labelfs=20, zoom=1, ticks=50,
               time_shift=0):

    if isinstance(inv, Inventory):
        attach_network_to_traces(st, inv)
        attach_coordinates_to_traces(st, inv, event)

        center = geometrical_center(inv)
        cstat = find_closest_station(inv, st, center['latitude'],
                                     center['longitude'])

        for i, trace in enumerate(st):
            if trace.stats.station not in [cstat]:
                continue
            else:
                sref = i
    else:
        sref = 0

    vespa = np.roll(data[0]*zoom, time_shift*int(st[0].stats.sampling_rate))
    taxis = data[1]
    urange = data[2]

    fig, ax = plt.subplots()
    try:
        refphase = st[0].stats.aligned
    except Exception:
        refphase = None
    if markphases:
        # RE = 6371.0
        try:
            origin = st[0].stats.origin
            depth = st[0].stats.depth

        except Exception:
            origin = event.origins[0]['time']
            depth = event.origins[0]['depth'] / 1000.

        m = TauPyModel('ak135')
        dist = st[sref].stats.distance
        arrival = m.get_travel_times(depth, dist, phase_list=markphases)

    # Labels of the plot.
    # Check if it is a relative plot to an aligned Phase.
    if refphase:
        try:
            m_tt = m.get_travel_times(depth, dist, refphase)
            p_ref = m_tt[0].ray_param_sec_degree
            ylabel = r'Relative $p$ in $\pm \frac{s}{deg}$  to %s arrival'
            ax.set_ylabel(ylabel % refphase, fontsize=fs)
            try:
                ax.set_title(r'Relative %ith root Vespagram' % (power),
                             fontsize=fs)
            except Exception:
                ax.set_title(r'Relative linear Vespagram', fontsize=fs)
        except Exception:
            p_ref = 0
            ylabel = r'Relative $p$ in $\pm \frac{s}{deg}$  to %s arrival'
            ax.set_ylabel(ylabel % refphase, fontsize=fs)
            try:
                ax.set_title(r'Relative %ith root Vespagram' % (power),
                             fontsize=fs)
            except Exception:
                ax.set_title(r'Relative linear Vespagram', fontsize=fs)
    else:
        p_ref = 0
        ax.set_ylabel(r'$p$ in $\frac{s}{deg}$', fontsize=fs)
        try:
            ax.set_title(r'%ith root Vespagram' % (power), fontsize=fs)
        except Exception:
            ax.set_title(r'Linear Vespagram', fontsize=fs)

    ax.set_xlabel(r'Time in s', fontsize=fs)

    # Do the contour plot of the Vespagram.
    if plot in ['contour', 'Contour']:
        if tw:
            tw = np.array(tw)
            cax = ax.imshow(vespa[tw.min():tw.max()], aspect='auto',
                            extent=(tw.min(), tw.max(),
                                    urange.min(), urange.max()),
                            origin='lower', cmap=cmap)
        else:
            cax = ax.imshow(vespa, aspect='auto',
                            extent=(taxis.min(), taxis.max(),
                                    urange.min(), urange.max()),
                            origin='lower', cmap=cmap)

        if markphases:
            for phase in arrival:
                t = phase.time
                phase_time = origin + t - st[sref].stats.starttime
                Phase_npt = int(phase_time / st[sref].stats.delta)
                tPhase = Phase_npt * st[sref].stats.delta
                name = phase.name
                sloPhase = phase.ray_param_sec_degree - p_ref
                if tPhase > taxis.max() or tPhase < taxis.min():
                    continue
                elif sloPhase > urange.max() or sloPhase < urange.min():
                    continue

                ax.autoscale(False)

                if marker in ['|']:
                    ax.plot((tPhase, tPhase),
                            (-labelfs/100. + sloPhase, labelfs/100.+sloPhase),
                            color=markerclr)
                else:
                    ax.plot(tPhase, sloPhase, marker)

                if name in [u'P^220P']:
                    ax.annotate('$P^{220}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

                elif name in [u'P^410P']:
                    ax.annotate('$P^{410}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

                elif name in [u'P^660P']:
                    ax.annotate('$P^{660}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase-0.6),
                                fontsize=labelfs, color=markerclr)

                else:
                    ax.annotate('$%s$' % name, xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

        cbar = fig.colorbar(cax, format='%.1f')
        cbar.set_clim(-1, 1)
        cbar.ax.tick_params(labelsize=fs)
        cbar.ax.set_ylabel('A', fontsize=fs)

    # Plot all the traces of the Vespagram.
    else:
        ax.set_ylim(urange[0] - 0.5, urange[urange.size - 1] + 0.5)
        ax.set_xticks(np.arange(taxis[0], taxis[taxis.size - 1], ticks))
        for i, trace in enumerate(vespa):
            ax.plot(taxis, trace + urange[i], color='black')

        if markphases:
            for phase in arrival:
                t = phase.time
                phase_time = origin + t - st[sref].stats.starttime
                Phase_npt = int(phase_time / st[sref].stats.delta)
                tPhase = Phase_npt * st[sref].stats.delta
                name = phase.name
                sloPhase = phase.ray_param_sec_degree - p_ref
                if tPhase > taxis.max() or tPhase < taxis.min():
                    continue
                elif sloPhase > urange.max() or sloPhase < urange.min():
                    continue

                if marker in ['|']:
                    ax.plot((tPhase, tPhase),
                            (-labelfs/100. + sloPhase, labelfs/100.+sloPhase),
                            color=markerclr)
                else:
                    ax.plot(tPhase, sloPhase, marker)

                if name in [u'P^220P']:
                    ax.annotate('$P^{220}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

                elif name in [u'P^410P']:
                    ax.annotate('$P^{410}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

                elif name in [u'P^660P']:
                    ax.annotate('$P^{660}P$', xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase-0.6),
                                fontsize=labelfs, color=markerclr)

                else:
                    ax.annotate('$%s$' % name, xy=(tPhase, sloPhase),
                                xytext=(tPhase + 1, sloPhase+0.2),
                                fontsize=labelfs, color=markerclr)

        if tw:
            tw = np.array(tw)
            pylab.xlim(tw.min(), tw.max())

    ax.tick_params(axis='both', which='major', labelsize=fs)

    if len(ax.xaxis.get_ticklabels()) > 5:
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)

    if savefig:
        fig.set_size_inches(10, 6)
        fig.savefig(savefig, dpi=dpi)
        plt.close("all")
    else:
        plt.ion()
        plt.draw()
        plt.show()
        plt.ioff()


def plot_vespa_stdout(data, st=None, inv=None, event=None, mp='P',
                      savefig=False, dpi=400, fs=25, power=4):

    if savefig:
        sf = savefig + "classic_no_captions.png"
        plot_vespa(data=data, st=st, inv=inv, event=event, markphases=False,
                   plot='classic', cmap='seismic', savefig=sf, dpi=dpi, fs=fs,
                   power=power)

        sf = savefig + "classic.png"
        plot_vespa(data=data, st=st, inv=inv, event=event,
                   markphases=[mp, 'P^410P', 'P^660P', 'PP'], plot='classic',
                   cmap='seismic', savefig=sf, dpi=dpi, fs=fs, power=power)

        sf = savefig + "contour_no_captions.png"
        plot_vespa(data=data, st=st, inv=inv, event=event, markphases=False,
                   plot='contour',
                   cmap='seismic', savefig=sf, dpi=dpi, fs=fs, power=power)

        sf = savefig + "contour.png"
        plot_vespa(data=data, st=st, inv=inv, event=event,
                   markphases=[mp, 'P^410P', 'P^660P', 'PP'], plot='contour',
                   cmap='seismic', savefig=sf, dpi=dpi, fs=fs, power=power)

    return


def mode_lines(ax, xrange, modes, overlap=True,
               label_height=None, label_width=None):
    xtrans = ax.get_xaxis_transform()
    mode_freqs = []
    mode_labels = []

    for mode in modes:
        if mode.freq > xrange[0] and mode.freq < xrange[-1]:
            name = r"$_{%s}%s_{%s}$" % (mode.n, mode.type, mode.l)
            mode_labels.append([mode.freq, name])
            mode_freqs.append(mode.freq)
    axvlines(mode_freqs, ax=ax, linestyle=':', linewidth=1,
             color='grey')
    if mode_labels:
        mode_labels = np.array(mode_labels, dtype=object)
        m_names = mode_labels.transpose()[1]

        if label_height is None or label_width is None:
            axtexts(mode_freqs, m_names,
                    y=1.03, ax=ax, transform=xtrans, rotation=45,
                    ha="center", va='center')
        else:
            m_ypos = np.ones(len(mode_freqs)) * ax.get_ylim()[1]
            txt_height = label_height*(ax.get_ylim()[1] - ax.get_ylim()[0])
            txt_width = label_width*(ax.get_xlim()[1] - ax.get_xlim()[0])
            text_positions = get_text_positions(mode_freqs, m_ypos, txt_width,
                                                txt_height)
            text_plotter(mode_freqs, m_ypos, text_positions, m_names, ax,
                         txt_width, txt_height)

    return


def plot_speclat(stream, faxis, tw, fw, zoom=1, segments=None, tr_ref=None,
                 highlight=None, modedata=None, ax=None,
                 **plot_kwargs):
    """
    if a segment file is given, tw is ignored
    """

    if highlight:
        if type(highlight) == str:
            highlight = [highlight]
    else:
        highlight = []

    if ax is None:
        fig, ax = plt.subplots()

    # Labels on axis
    ax.set_ylabel('Latitude (deg)')
    ax.set_xlabel('frequency (mHz)')

    """
    Create one polygon patch here for each segment-collection per mode
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    patches []
    polygon = Polygon(np.array([[fw1, lat1], [fw1, lat2],
                                [fw2, lat1], [fw2, lat2]]), True)
    patches.append(polygon)
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

    # Set up colorscale
    colors = 100*np.random.rand(len(patches))
    p.set_array(np.array(colors))

    ax.add_collection(p)
    """
    # Loop for max amplitude in stream:
    if zoom == 'normalize each':
        max_amp = 0.0
        for tr in stream:
            if segments:
                seg = segments.select(station=tr.stats.station)[0]
                tw = [seg.tw1, seg.tw2]

                if fw is None:
                    fw = [seg.fw1, seg.fw2]
            f, Fxx, startlabel, endlabel = taper_FT(tr, tw, fw)
            amp = abs(Fxx)
            if amp[startlabel:endlabel+1].max() > max_amp:
                max_amp = amp[startlabel:endlabel+1].max()
        zoom = 1./max_amp * 20

    # Loop for spectra
    amps = None
    freqs = None

    amps_hl = None
    freqs_hl = None

    lat_label = AttribDict()

    for i, tr in enumerate(stream):
        if segments:
            try:
                seg = segments.select(station=tr.stats.station)[0]
                tw = [seg.tw1, seg.tw2]
            except IndexError:
                continue

        f, Fxx, startlabel, endlabel = taper_FT(tr, tw, fw)
        amp = abs(Fxx)

        lat = tr.stats.ah.station.latitude

        freq = np.append(f[startlabel:endlabel+1], np.nan)
        amp_tr = amp[startlabel:endlabel+1]
        if type(zoom) == float or type(float) == int:
            amp = np.append(amp_tr*float(zoom)/amp_tr.max() + lat, np.nan)
        else:
            amp = np.append(zoom*amp_tr + lat, np.nan)

        #
        if tr.stats.station in highlight or tr == tr_ref:
            if freqs_hl is None:
                freqs_hl = freq
                amps_hl = amp
            else:
                freqs_hl = np.vstack((freqs, freq))
                amps_hl = np.vstack((amps, amp))
        else:
            if freqs is None:
                freqs = freq
                amps = amp
            else:
                freqs = np.vstack((freqs, freq))
                amps = np.vstack((amps, amp))
        lat_label[tr.stats.station] = lat

    if freqs is not None:
        freqs = freqs.flatten()
        amps = amps.flatten()

        # Plot spectra
        ax.plot(freqs, amps, **plot_kwargs)

    if freqs_hl is not None:
        ax.plot(freqs_hl, amps_hl, color='red')

    # Loop for modes
    if modedata is not None:
        mode_lines(ax, faxis, modedata)

    # Loop for station names
    ytrans = ax.get_yaxis_transform()
    for station, lat in lat_label.items():
        ax.annotate("%s" % station, xy=(0.05, lat), xycoords=ytrans)

    plt.show()
    return


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.

    H. Karaoglu
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def get_iter_colormap(input_list, cmap, random_values=False):
    if cmap == 'GreensBlues' and len(input_list) > 3:
        Greens = cm.get_cmap('Greens', 256)
        Blues = cm.get_cmap('Blues', 256)
        green = Greens(np.linspace(0.25, 0.75, len(input_list)/2))
        blue = Blues(np.linspace(
                        0.25, 0.75, len(input_list)-len(input_list)/2
                        ))
        blue = np.flip(blue, 0)
        GreenBlue = np.concatenate((blue, green))
        colormap = iter(GreenBlue)
    elif cmap == 'GreensBluesBlack' and len(input_list) > 3:
        Greens = cm.get_cmap('Greens', 256)
        Blues = cm.get_cmap('Blues', 256)
        Blacks = cm.get_cmap('Greys', 256)
        green = Greens(np.linspace(0.25, 0.75, len(input_list)/2))
        blue = Blues(np.linspace(0.25, 0.75,
                                 len(input_list)-len(input_list)/2-1))
        blue = np.flip(blue, 0)
        black = Blacks(np.linspace(0.6, 0.6, 1))  # (np.linspace(1, 1, 1))
        GreenBlueBlack = np.concatenate((blue, green, black))
        colormap = iter(GreenBlueBlack)
    elif cmap == 'BlackGreysRed' and len(input_list) > 1:
        Reds = cm.get_cmap('Reds', 256)
        Grays = cm.get_cmap('binary', 256)
        Blacks = cm.get_cmap('Greys', 256)
        red = Reds(np.linspace(0.6, 0.6, 1))
        black = Blacks(np.linspace(1, 1, 1))
        gray = Grays(np.linspace(0.6, 0.8, len(input_list)-2))
        BlackGreysRed = np.concatenate((black, gray, red))
        colormap = iter(BlackGreysRed)
    elif cmap == 'BlackRedGreys' and len(input_list) > 1:
        Reds = cm.get_cmap('Reds', 256)
        Grays = cm.get_cmap('binary', 256)
        Blacks = cm.get_cmap('Greys', 256)
        red = Reds(np.linspace(0.6, 0.6, 1))
        black = Blacks(np.linspace(1, 1, 1))
        gray = Grays(np.linspace(0.6, 0.8, len(input_list)-2))
        BlackRedGreys = np.concatenate((black, red, gray))
        colormap = iter(BlackRedGreys)
    elif cmap == 'BlueBlackGreysRed' and len(input_list) > 1:
        Reds = cm.get_cmap('Reds', 256)
        Grays = cm.get_cmap('binary', 256)
        Blacks = cm.get_cmap('Greys', 256)
        Blues = cm.get_cmap('winter', 256)
        red = Reds(np.linspace(0.6, 0.6, 1))
        black = Blacks(np.linspace(1, 1, 1))
        blue = Blues(np.linspace(0, 0, 1))
        gray = Grays(np.linspace(0.6, 0.8, len(input_list)-3))
        BlueBlackGreysRed = np.concatenate((blue,black,gray,red))
        colormap = iter(BlueBlackGreysRed)
    elif cmap == 'Grays' and len(input_list) > 1:
        gray = cm.get_cmap('binary', 256)
        Grays = gray(np.linspace(0.25, 0.75, len(input_list)))
        colormap = iter(Grays)
    elif cmap == 'tab10' and random_values is False:
        colormap = iter(getattr(cm, cmap)(
                        np.linspace(0, 1, 10))
                        )
    elif cmap == 'Set1' and random_values is False:
        colormap = iter(getattr(cm, cmap)(
                        np.linspace(0, 1, 9))
                        )

    elif random_values != False and random_values != 'center':
        colormap = iter(getattr(cm, cmap)(
                        np.random.rand(len(input_list)))
                        )
    elif random_values == 'center':
        colormap = iter(getattr(cm, cmap)(
                        np.linspace(0.1, 0.9, len(input_list)))
                        )
    elif len(input_list) > 3:
        colormap = iter(getattr(cm, cmap)(
                        np.linspace(0, 1, len(input_list)))
                        )
    else:
        colormap = iter(['blue', 'red', 'green'])
    return colormap


def get_iter_linestyle(exclude=None):
    lines = ["-", ":", "--", "-."]
    if exclude is not None:
        lines.remove(exclude)
    linecycler = cycle(lines)
    return linecycler


def get_text_positions(x_data, y_data, txt_width, txt_height):
    text_positions = y_data.copy()

    # ylimit = text_positions.max()
    # ydiff = 0

    for index, (y, x) in enumerate(zip(y_data, x_data)):
        local_tp = [i for i in zip(y_data, x_data) if i[0] > (y - txt_height)
                    and (abs(i[1] - x) < txt_width * 2) and i != (y, x)]
        if local_tp:
            sorted_ltp = sorted(local_tp)
            if abs(sorted_ltp[0][0] - y) < txt_height:  # True == collision
                differ = np.diff(sorted_ltp, axis=0)
                y_data[index] = sorted_ltp[-1][0] + txt_height
                text_positions[index] = sorted_ltp[-1][0] + txt_height
                for k, (j, m) in enumerate(differ):
                    # j is the vertical distance between words
                    if j > txt_height * 2:  # if True then room to fit a word
                        y_data[index] = sorted_ltp[k][0] + txt_height
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    #     if text_positions[index] > ylimit:
    #         _ydiff = text_positions[index] - ylimit
    #         if _ydiff > ydiff:
    #             ydiff = _ydiff
    # if ydiff > 0:
    #     text_positions = text_positions - ydiff
    return text_positions


def text_plotter(x_data, y_data, text_positions, text_label, axis, txt_width,
                 txt_height, fontsize=10):
    """
    #random test data:
    x_data = random_sample(100)
    y_data = random_integers(10,50,(100))

    #GOOD PLOT:
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.bar(x_data, y_data,width=0.00001)
    #set the bbox for the text. Increase txt_width for wider text.
    txt_height = 0.04*(plt.ylim()[1] - plt.ylim()[0])
    txt_width = 0.02*(plt.xlim()[1] - plt.xlim()[0])
    #Get the corrected text positions, then write the text.
    text_positions = get_text_positions(x_data, y_data, txt_width, txt_height)
    text_plotter(x_data, y_data, text_positions, ax2, txt_width, txt_height)

    plt.ylim(0,max(text_positions)+2*txt_height)
    plt.xlim(-0.1,1.1)

    #BAD PLOT:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(x_data, y_data, width=0.0001)
    #write the text:
    for x,y in zip(x_data, y_data):
        ax.text(x - txt_width, 1.01*y, '%d'%int(y),rotation=0)
    plt.ylim(0,max(text_positions)+2*txt_height)
    plt.xlim(-0.1,1.1)

    plt.show()
    """
    ylim = axis.get_ylim()
    for x, y, t, text in zip(x_data, y_data, text_positions, text_label):
        axis.text(x, 1.01 * t, '%s' % text, rotation=0,
                  color='black', ha="center", va='center', fontsize=fontsize)
        if y != t:
            axis.plot([x, x], [t, y], color='grey', linestyle=':',
                      linewidth=1)
        #     axis.arrow(x, t, 0, y-t, color='blue', alpha=0.3,
        #                width=txt_width*0.2, head_width=txt_width,
        #                head_length=txt_height*0.5,
        #                zorder=0, length_includes_head=True)

    axis.set_ylim(ylim[0], max(text_positions) * 1.1)


def plot_Mmatrix(matrix_file, omega_file, modes_dir, matrix_type='bin',
                 title=None, show=True, Qoption=False, savefig=False, **kwargs):
    """
    Prints M matrix, real and imaginary parts together with singlets (f vs Q)
    matrix_type: bin or ascii, if ascii first line is matrix size
    Qoption: default we plot Q, optional to plot 1000/Q
    """
    if matrix_type == 'bin':
        N, M_re, M_im = bin2ascii_matrix(matrix_file)  # N: matrix size
        M_re = M_re/1e-6  # converting to muHz
        M_im = M_im/1e-6  # converting to muHz
    else:
        N = int(open(matrix_file).readline().rstrip())  # matrix size
        M = np.loadtxt(matrix_file, skiprows=1)
        M = [ii for ii in zip(*M)]
        M_re = np.array(M[0]).reshape((N, -1))/1e-6
        M_im = np.array(M[1]).reshape((N, -1))/1e-6

    omega = read_omega_dat(omega_file)
    omega = [ii for ii in zip(*omega)]
    modes, modes_cc = read_modes_in(modes_dir)
    modes_all = read_modes()
    n = int(modes[0])  # number of modes in file

    # colormap
    vminp = min(np.min(M_re), np.min(M_im))
    vmaxp = np.max(M_im)

    if 'vmin' in kwargs:
        vminp = kwargs['vmin']
    if 'vmax' in kwargs:
        vmaxp = kwargs['vmax']

    if vminp < 0 and vmaxp > 0:
        mpnt = 0.
    else:
        mpnt = 1 - vmaxp/(vmaxp + abs(vminp))

    colors = [(0, 0, 1), (0.9, 0.9, 0), (1, 0, 0)]  # B -> Y -> R
    cmap_name = 'sf_cmap'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=7)
    cmap = shiftedColorMap(cmap, start=0,
                           midpoint=(mpnt-vminp)/(vmaxp-vminp),
                           stop=1.0, name='shiftedcmap')

    fig = plt.figure(figsize=(15, 7))
    ax = []

    # f vs Q
    ax.append(fig.add_subplot(1, 3, 1))
    if Qoption:
        omega[1] = [1000./x for x in omega[1]]
    ax[0].scatter(omega[0], omega[1], color='r', marker='D', s=10)
    for mode in modes[1:n+1]:
        mode = modes_all.select(name=''.join(mode.split()[:3]))[0]
        label = '${}_{%s}%s_{%s}$' % (mode.n, mode.type, mode.l)
        ax[0].text(mode.freq, mode.Q, label, fontsize=16)
    if 'xlim' in kwargs:
        ax[0].set_xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        ax[0].set_ylim(kwargs['ylim'])
    aspect = (ax[0].get_xlim()[1] - ax[0].get_xlim()[0]) /   \
             (ax[0].get_ylim()[1] - ax[0].get_ylim()[0])
    ax[0].set_aspect(adjustable='box', aspect=aspect)
    ax[0].set_xlabel("f (mHz)", fontsize=16)
    if Qoption:
        ax[0].set_ylabel("1000/Q", fontsize=16)
    else:
        ax[0].set_ylabel("Q", fontsize=16)

    # Re[M]
    ax.append(fig.add_subplot(1, 3, 2))
    s = ax[1].imshow(M_re, cmap=cmap, aspect=1, vmin=vminp, vmax=vmaxp)
    mf = ticker.ScalarFormatter(useMathText=True)  # colorbar
    mf.set_powerlimits((-2, 2))  # for exponents greater than +-2
    cb = fig.colorbar(s, orientation='horizontal', pad=0.01, format=mf)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.set_xlabel('Re[$M$] ($\mu$Hz)', fontsize=16)
    ax[1].tick_params(axis='x', which='both', top='on', bottom='off',
                      labelbottom='off', labeltop='on')

    lsum = 0
    l = []
    m = []
    d = 0.5  # d: size of point in imshow
    for mode in modes[1:n+1]:  # creating lines in matrix
        mode = modes_all.select(name=''.join(mode.split()[:3]))[0]
        lsum = lsum + (2*mode.l + 1)
        m.extend(range(-mode.l, mode.l+1))  # tick positions and singlet names
        ax[1].plot([lsum-d, lsum-d], [0-d, N-d], lw=1, c="w")  # vertical lines
        ax[1].plot([0-d, N-d], [lsum-d, lsum-d], lw=1, c="w")  # hor. lines
        l.append(mode.l)

    if max(l) <= 12:  # printing singlet names
        ax[1].set_xticks(range(len(m)))
        ax[1].set_xticklabels(m, fontsize=6)
        ax[1].set_yticks(range(len(m)))
        ax[1].set_yticklabels(m, fontsize=6)
        ax[1].set_xlabel('$m$', fontsize=16)
        ax[1].set_ylabel('$m^{\prime}$', rotation='vertical', fontsize=16)
        ax[1].xaxis.set_label_position('top')
        ax[1].set_xlim(-d, N-1+d)
        ax[1].set_ylim(N-1+d, -d)

    # Im[M]
    ax.append(fig.add_subplot(1, 3, 3))
    s = ax[2].imshow(M_im, cmap=cmap, aspect=1, vmin=vminp, vmax=vmaxp)
    mf = ticker.ScalarFormatter(useMathText=True)  # colorbar
    mf.set_powerlimits((-2, 2))  # for exponents greater than +-2
    cb = fig.colorbar(s, orientation='horizontal', pad=0.01, format=mf)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.set_xlabel('Im[$M$] ($\mu$Hz)', fontsize=16)
    ax[2].tick_params(axis='x', which='both', top='on', bottom='off',
                      labelbottom='off', labeltop='on')

    lsum = 0
    l = []
    m = []
    d = 0.5  # d: size of point in imshow
    for mode in modes[1:n+1]:  # creating lines in matrix
        mode = modes_all.select(name=''.join(mode.split()[:3]))[0]
        label = '${}_{%s}%s_{%s}$' % (mode.n, mode.type, mode.l)
        lsum = lsum + (2*mode.l + 1)
        m.extend(range(-mode.l, mode.l+1))  # tick positions and singlet names
        ax[2].plot([lsum-d, lsum-d], [0-d, N-d], lw=1, c="w")  # vertical lines
        ax[2].plot([0-d, N-d], [lsum-d, lsum-d], lw=1, c="w")  # hor. lines
        ax[2].annotate(label, xy=(2, 1), xytext=(N-1+d, lsum-mode.l-d),
                       va='center', ha='left', fontsize=16)
        l.append(mode.l)

    if max(l) <= 12:  # printing singlet names
        ax[2].set_xticks(range(len(m)))
        ax[2].set_xticklabels(m, fontsize=6)
        ax[2].set_yticks(range(len(m)))
        ax[2].set_yticklabels(m, fontsize=6)
        ax[2].set_xlabel('$m$', fontsize=16)
        ax[2].set_ylabel('$m^{\prime}$', rotation='vertical', fontsize=16)
        ax[2].xaxis.set_label_position('top')
        ax[2].set_xlim(-d, N-1+d)
        ax[2].set_ylim(N-1+d, -d)

    plt.show()
    if savefig:
        plt.tight_layout()
        fig.savefig('matrix.png', orientation='landscape', dpi=400,
                    bbox_inches='tight', pad_inches=0.025)
    return


def freq_over_l(mtype='T', lmin=2, lmax=40, nmin=0, nmax=None,
                fmin=None, fmax=None,
                boxes=None, plot_boxes_as_points=False, boxmarker="o",
                seg_scaling=0.5, verbose=False,
                highlight=[None], deltax_text=0.25, deltay_text=0,
                show=False, fig_size=(7, 12), font_size=22,
                tick_width=1, border_width=1, loc='best',
                xlim=None, ylim=None, plot_text=True,
                **pltargs):
    """
    seg_count: dict, with modename as keys and number of segments as values
    highlight: dict, containing l or n as key and color as value

    m = {
     '00t05': 100,
     '00t06': 100,
     '00t07': 100,
     }

    For presentation
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    ax.xaxis.label.set_fontsize(20)
    ax.yaxis.label.set_fontsize(20)
    """
    if mtype not in ["T", "S", "TZ", "cc"]:
        print('mtype must be "T", "S" or "cc"')
        return

    def plot_branch(branch, freqs, n, ls, color, marker, label):
        if n in highlight:
            ax.plot(branch, freqs, '-o', color=highlight[n],
                    label=label)
            if plot_text is True:
                ax.text(lmin - deltay_text, freqs[0] + deltay_text, "n=%i" % n)
        else:
            ax.plot(branch, freqs, linestyle=ls, marker=marker, color=color,
                    label=label)
            if plot_text is True:
                ax.text(lmin - deltax_text, freqs[0] + deltay_text, "n=%i" % n)

    def plot_type(modes, ls, color, marker, label, seg_box):
        n = 0
        branch = []
        freqs = []
        seg = []
        seg_freq = []
        seg_size = []
        mname = []
        for m in modes:
            if m.n < nmin:
                continue

            if nmax is not None:
                if m.n > nmax:
                    continue

            if m.l > lmax or m.l < lmin:
                continue

            if fmin is not None:
                if m.freq < fmin:
                    continue

            if fmax is not None:
                if m.freq > fmax:
                    continue

            if m.n != n:
                if len(branch) != 0:
                    plot_branch(branch, freqs, n, ls, color, marker, label)
                    if verbose is True:
                        print(mname)
                n = m.n
                branch = []
                freqs = []
                mname = []

            mname += [m.name]
            branch += [m.l]
            freqs += [m.freq]

            if boxes is not None:
                if any(m.folder_name() in s for s in boxes.keys()):
                    seg += [m.l]
                    seg_freq += [m.freq]
                    keys = list(filter(lambda x: m.folder_name() in x, boxes))
                    seg_size += [boxes[keys[0]]]
                    if seg_box is not None:
                        for k in keys:
                            if k not in seg_box:
                                seg_box[k] = [[m.l, m.freq]]
                            else:
                                x = seg_box[k]
                                x.append([m.l, m.freq])
                                seg_box[k] = x
        return branch, freqs, n, mname, seg_size, seg, seg_freq, seg_box

    all_modes = read_modes()

    if mtype in ['cc', 'TZ']:
        types = ["S", "T"]
        cc_box = {}
        # color = {'S': '#DC143C', 'T': '#20B2AA'}
        color = {'S': '#616161', 'T': '#515151'}
        ls = {'S': '-', 'T': '--'}
        marker = {'S': 'o', 'T': '^'}
    else:
        cc_box = None
        color = {mtype.upper(): 'black'}
        marker = {mtype.upper(): 'o'}
        ls = {mtype.upper(): '-'}
        types = [mtype.upper()]

    label = {'S': 'Spheroidal Modes', 'T': 'Toroidal Modes'}
    fig, ax = plt.subplots()
    for i, mt in enumerate(types):
        modes = all_modes.select(mtype=mt)
        modes.sort(['n', 'l'])

        out = plot_type(modes, ls[mt], color[mt], marker[mt],
                        label=None, seg_box=cc_box)
        branch, freqs, n, mname, seg_size, seg, seg_freq, cc_box = out[:]
        if len(branch) != 0:
            plot_branch(branch, freqs, n, ls[mt], color[mt], marker[mt],
                        label=label[mt])
            if verbose is True:
                print(mname)

        if boxes is not None:
            if verbose:
                print('Angular orders l:')
                print(seg)
            seg_size = [x*seg_scaling for x in seg_size]
            if i == 0:
                if plot_boxes_as_points is False:
                    ax.scatter(seg, seg_freq, s=seg_size, marker='s',
                               facecolors='none', edgecolors='r', zorder=100,
                               label='Measured Modes')
                else:
                    ax.plot(seg, seg_freq, linestyle=None, linewidth=0,
                            marker=marker[mt], #boxmarker,
                            color='red',
                            label='Measured Modes')
            else:
                if plot_boxes_as_points is False:
                    ax.scatter(seg, seg_freq, s=seg_size, marker='s',
                               facecolors='none', edgecolors='r', zorder=100)
                else:
                    ax.plot(seg, seg_freq, linestyle=None, linewidth=0,
                            marker=marker[mt], #boxmarker,
                            color='red')

    if cc_box is not None:
        for i, val in enumerate(cc_box.values()):
            # x = [val[0][0]-0.1, val[0][0]-0.1,
            #      val[1][0]+0.1, val[1][0]+0.1,
            #      val[0][0]-0.1]
            # y = [val[0][1]-0.05, val[0][1]+0.05,
            #      val[1][1]+0.05, val[1][1]-0.05,
            #      val[0][1]-0.05]
            x, y = np.array(val).transpose()
            if i == 0:
                ax.plot(x, y, color='r', zorder=1, label="Cross-coupling")
            else:
                ax.plot(x, y, color='r', zorder=1)

    fig.set_size_inches(fig_size)
    ax.set_ylabel('frequency (mHz)')
    ax.set_xlabel('Angular order l')

    if types == ["T"]:
        title = "Toroidal Modes"
    elif types == ["S"]:
        title = "Spheroidal Modes"
    elif types == ["S", "T"]:
        title = "Cross-Coupling"

    ax.legend(loc=loc)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_tick_params(which='both', width=tick_width)
    ax.yaxis.set_tick_params(which='both', width=tick_width)
    ax.set_title(title)
    mpl.rcParams.update({'font.size': font_size})

    [i.set_linewidth(border_width) for i in ax.spines.values()]

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if show is True:
        plt.show()
    return fig, ax
