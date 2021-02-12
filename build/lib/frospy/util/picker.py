from __future__ import absolute_import, print_function

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mpld
import matplotlib.cbook as cbook
import matplotlib.gridspec as gridspec

from matplotlib.lines import Line2D
from matplotlib import path as mplPath
from matplotlib.widgets import Button, CheckButtons

import numpy as np
try:
    import scipy.spatial as spatial
    SCIPY = True
except ImportError:
    SCPY = False

import datetime as dt
import sys

from frospy.plot.nmplt import format_exponent
from frospy.util.base import pairwise, mask_data, fourier_transform


def fmt(x, y):
    return 'x: {x:0.2f}\ny: {y:0.2f}'.format(x=x, y=y)


class DataCursor(object):
    # http://stackoverflow.com/a/4674445/190597
    """A simple data cursor widget that displays the x,y location of a
    matplotlib artist when it is selected.
    Example:

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    scat = ax.plot(x, y)
    DataCursor(scat, x, y)
    plt.show()


    """
    def __init__(self, artists, x=[], y=[], tolerance=5, offsets=(-20, 20),
                 formatter=fmt, display_all=False):
        """Create the data cursor and connect it to the relevant figure.
        "artists" is the matplotlib artist or sequence of artists that will be
            selected.
        "tolerance" is the radius (in points) that the mouse click must be
            within to select the artist.
        "offsets" is a tuple of (x,y) offsets in points from the selected
            point to the displayed annotation box
        "formatter" is a callback function which takes 2 numeric arguments and
            returns a string
        "display_all" controls whether more than one annotation box will
            be shown if there are multiple axes.  Only one will be shown
            per-axis, regardless.
        """
        self._points = np.column_stack((x, y))
        self.formatter = formatter
        self.offsets = offsets
        self.display_all = display_all
        if not cbook.iterable(artists):
            artists = [artists]
        self.artists = artists
        self.axes = tuple(set(art.axes for art in self.artists))
        self.figures = tuple(set(ax.figure for ax in self.axes))

        self.annotations = {}
        for ax in self.axes:
            self.annotations[ax] = self.annotate(ax)

        for artist in self.artists:
            artist.set_picker(tolerance)
        for fig in self.figures:
            fig.canvas.mpl_connect('pick_event', self)

    def annotate(self, ax):
        """Draws and hides the annotation box for the given axis "ax"."""
        annotation = ax.annotate(self.formatter, xy=(0, 0), ha='right',
                                 xytext=self.offsets,
                                 textcoords='offset points', va='bottom',
                                 bbox=dict(boxstyle='round,pad=0.5',
                                           fc='yellow', alpha=0.5),
                                 arrowprops=dict(arrowstyle='->',
                                                 connectionstyle='arc3,rad=0')
                                 )
        annotation.set_visible(False)
        return annotation

    def snap(self, x, y):
        """Return the value in self._points closest to (x, y).
        """
        idx = np.nanargmin(((self._points - (x, y))**2).sum(axis=-1))
        return self._points[idx]

    def __call__(self, event):
        """Intended to be called through "mpl_connect"."""
        # Rather than trying to interpolate, just display the clicked coords
        # This will only be called if it's within "tolerance", anyway.
        x, y = event.mouseevent.xdata, event.mouseevent.ydata
        annotation = self.annotations[event.artist.axes]
        if x is not None:
            if not self.display_all:
                # Hide any other annotation boxes...
                for ann in self.annotations.values():
                    ann.set_visible(False)
            # Update the annotation in the current axis..
            x, y = self.snap(x, y)
            annotation.xy = x, y
            annotation.set_text(self.formatter(x, y))
            annotation.set_visible(True)
            event.canvas.draw()


class FollowDotCursor(object):
    """
    Display the x,y location of the nearest data point.

    Example:

    x=[1,2,3,4,5]
    y=[6,7,8,9,10]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x, y)
    cursor = FollowDotCursor(ax, x, y)
    plt.show()

    """
    def __init__(self, ax, x, y, tolerance=5, formatter=fmt,
                 offsets=(-20, 20)):
        try:
            x = np.asarray(x, dtype='float')
        except (TypeError, ValueError):
            x = np.asarray(mpld.date2num(x), dtype='float')
        y = np.asarray(y, dtype='float')
        self._points = np.column_stack((x, y))
        self.offsets = offsets
        self.scale = x.ptp()
        self.scale = y.ptp() / self.scale if self.scale else 1
        self.tree = spatial.cKDTree(self.scaled(self._points))
        self.formatter = formatter
        self.tolerance = tolerance
        self.ax = ax
        self.fig = ax.figure
        self.ax.xaxis.set_label_position('top')
        self.dot = ax.scatter(
            [x.min()], [y.min()], s=130, color='green', alpha=0.7)
        self.annotation = self.setup_annotation()
        plt.connect('motion_notify_event', self)

    def scaled(self, points):
        points = np.asarray(points)
        return points * (self.scale, 1)

    def __call__(self, event):
        ax = self.ax
        # event.inaxes is always the current axis. If you use twinx, ax could be
        # a different axis.
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
        elif event.inaxes is None:
            return
        else:
            inv = ax.transData.inverted()
            x, y = inv.transform([(event.x, event.y)]).ravel()
        annotation = self.annotation
        x, y = self.snap(x, y)
        annotation.xy = x, y
        annotation.set_text(self.formatter(x, y))
        self.dot.set_offsets((x, y))
        bbox = ax.viewLim
        event.canvas.draw()

    def setup_annotation(self):
        """Draw and hide the annotation box."""
        annotation = self.ax.annotate(
            '', xy=(0, 0), ha = 'right',
            xytext = self.offsets, textcoords = 'offset points', va = 'bottom',
            bbox = dict(
                boxstyle='round,pad=0.5', fc='yellow', alpha=0.75),
            arrowprops = dict(
                arrowstyle='->', connectionstyle='arc3,rad=0'))
        return annotation

    def snap(self, x, y):
        """Return the value in self.tree closest to x, y."""
        dist, idx = self.tree.query(self.scaled((x, y)), k=1, p=1)
        try:
            return self._points[idx]
        except IndexError:
            # IndexError: index out of bounds
            return self._points[0]

def get_polygon(data, no_of_vert=4, xlabel=None, xticks=None, ylabel=None, yticks=None, fs=25):
    """
    Interactive function to pick a polygon out of a figure and receive the vertices of it.
    :param data:
    :type:

    :param no_of_vert: number of vertices, default 4,
    :type no_of_vert: int
    """
    from frospy.util.polygon_interactor import PolygonInteractor
    from matplotlib.patches import Polygon

    no_of_vert = int(no_of_vert)
    # Define shape of polygon.
    try:
        x, y = xticks.max(), yticks.max()
        xmin= -x/10.
        xmax= x/10.
        ymin= y - 3.*y/2.
        ymax= y - 3.*y/4.

    except AttributeError:
        y,x = data.shape
        xmin= -x/10.
        xmax= x/10.
        ymin= y - 3.*y/2.
        ymax= y - 3.*y/4.

    xs = []
    for i in range(no_of_vert):
        if i >= no_of_vert/2:
            xs.append(xmax)
        else:
            xs.append(xmin)

    ys = np.linspace(ymin, ymax, no_of_vert/2)
    ys = np.append(ys,ys[::-1]).tolist()

    poly = Polygon(list(zip(xs, ys)), animated=True, closed=False, fill=False)

    # Add polygon to figure.
    fig, ax = plt.subplots()
    ax.add_patch(poly)
    p = PolygonInteractor(ax, poly)
    plt.title("Pick polygon, close figure to save vertices")
    plt.xlabel(xlabel, fontsize=fs)
    plt.ylabel(ylabel, fontsize=fs)

    try:
        im = ax.imshow(abs(data), aspect='auto', extent=(xticks.min(), xticks.max(), 0, yticks.max()), interpolation='none')
    except AttributeError:
        im = ax.imshow(abs(data), aspect='auto', interpolation='none')

    cbar = fig.colorbar(im)
    cbar.ax.tick_params(labelsize=fs)
    cbar.ax.set_ylabel('R', fontsize=fs)
    mpl.rcParams['xtick.labelsize'] = fs
    mpl.rcParams['ytick.labelsize'] = fs
    ax.tick_params(axis='both', which='major', labelsize=fs)


    plt.show()
    print("Calculate area inside chosen polygon\n")
    try:
        vertices = (poly.get_path().vertices)
        vert_tmp = []
        xticks.sort()
        yticks.sort()
        for fkvertex in vertices:
            vert_tmp.append([np.abs(xticks-fkvertex[0]).argmin(), np.abs(yticks[::-1]-fkvertex[1]).argmin()])
        vertices = np.array(vert_tmp)

    except AttributeError:
        vertices = (poly.get_path().vertices).astype('int')

    indicies = convert_polygon_to_flat_index(data, vertices)
    return indicies


def convert_polygon_to_flat_index(data, vertices):
    """
    Converts points insde of a polygon defined by its vertices, taken of an
    imshow plot of data,to flat-indicies. Does NOT include the border of the
    polygon.

    :param data: speaks for itself
    :type data: numpy.ndarray

    :param vertices: also...
    :type vertices: numpy.ndarray

    """

    # check if points are inside polygon. Be careful with the indicies, np and
    # mpl handle them exactly opposed.
    polygon = mplPath.Path(vertices)
    arr = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if polygon.contains_point([j, i]):
                arr.append([j, i])
    arr = map(list, zip(*arr))

    flat_index = np.ravel_multi_index(arr, data.conj().transpose().shape)
    flat_index = flat_index.astype('int').tolist()

    return(flat_index)


def pick_data(x, y, ax=None, xlabel=None, ylabel=None, title=None,
              cursor=False, singlepick=False):

    xlim = None
    if not ax:
        fig, ax1 = plt.subplots(1, 1)
        ax1.set_title(title)
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)

    else:
        fig = plt.gcf()
        ax1 = ax
        ax1.clear()
        xlim = ax1.get_xlim()

    ax1.clear()
    ax1.plot(x, y, picker=5, color='black')  # 5 points tolerance
    ax1 = format_exponent(ax1)
    if xlim:
        ax1.set_xlim(xlim)

    if cursor and SCIPY is True:
        FollowDotCursor(ax1, x, y)
    # global PickByHand
    PickByHand = []

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            pick_tmp = zip(np.take(xdata, ind), np.take(ydata, ind))
            amax = np.array(pick_tmp).transpose()[1].argmax()
            pick = pick_tmp[amax]

            ax1.axvline(x=pick[0], linestyle='-',
                        linewidth=1, color='blue')

            PickByHand.append(pick)
            fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick1)

    if type(x[0]) == dt.datetime:
        plt.xticks(rotation=25)
        ax = plt.gca()
        xfmt = mpld.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.xaxis.set_major_formatter(xfmt)

    plt.show()
    return PickByHand


def pick_segment_window(x, y, ysyn=None, ax=None, xlabel=None, ylabel=None,
                        ytick_label=None, yticks=None,
                        xtick_label=None, xticks=None,
                        title=None, cursor=False, pairwisepick=True,
                        modes=None):

    class Index(object):
        ind = 0

        def delete(self, event):
            if len(pick_lines) != 0:
                PickByHand.pop()
                pick_lines[-1].remove()
                pick_lines.pop()
                plt.draw()
            else:
                return

    xlim = None
    if not ax:
        fig, ax1 = plt.subplots(1, 1)
        ax1.set_title(title)
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)

    else:
        fig = plt.gcf()
        ax1 = ax
        ax1.clear()
        xlim = ax1.get_xlim()

    l, = ax1.plot(x, y, picker=5, color='black')  # 5 points tolerance
    if ysyn is not None:
        _l, = ax1.plot(x, ysyn, color='blue', linestyle='dashed')
    ax1 = format_exponent(ax1)
    fig.set_size_inches(14, 10)
    if modes is not None:
        trans = ax1.get_xaxis_transform()
        for mode in modes:
            if mode.freq > x[0] and mode.freq < x[-1]:
                ax1.axvline(x=mode.freq, linestyle=':',
                            linewidth=1, color='grey')
                ax1.text(mode.freq, 1.03, r"$%s$" % mode.name,
                         transform=trans, rotation=45, ha="center",
                         va='center')
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)

    if xtick_label is not None:
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xtick_label)

    if xlim:
        ax1.set_xlim(xlim)

    if cursor and SCIPY is True:
        FollowDotCursor(ax1, x, y)
    global PickByHand
    PickByHand = []
    pick_lines = []

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            pick_tmp = zip(np.take(xdata, ind), np.take(ydata, ind))
            amax = np.array(pick_tmp).transpose()[1].argmax()
            pick = pick_tmp[amax]
            p = ax1.axvline(x=pick[0], linestyle='-',
                            linewidth=1, color='blue')
            pick_lines.append(p)
            plt.draw()

            PickByHand.append(pick[0])

    fig.canvas.mpl_connect('pick_event', onpick1)

    callback = Index()

    axdelete = plt.axes([0.91, 0.35, 0.08, 0.075])

    bdelete = Button(axdelete, 'Remove Pick')
    bdelete.on_clicked(callback.delete)

    plt.show()
    input('Press enter, when picked...')
    if pairwisepick:
        if len(PickByHand) % 2 != 0:
            msg = 'Number of picks must be multiple of 2, no window selected'
            print(msg)
            return []

    PickByHand.sort()
    return PickByHand


def pick_seismo_window(trace, trace_org, windows, xlabel, ylabel, title,
                       ylim=None, autoscale=False, cursor=False,
                       yaxis='dates'):

    def choose_window(data, dates, i, windows):
        win_trace = data[windows[i]:windows[i+1]]
        win_dates = dates[windows[i]:windows[i+1]]
        return win_dates, win_trace

    class Index(object):
        ind = 0

        def next(self, event):
            if len(PickByHand) % 2 == 0:
                self.ind += 1
                noi = len(windows)-1
                i = self.ind % noi
                if i == 0:
                    msg = 'First timewindow of station: '
                    ax1.set_title(msg + title)
                else:
                    ax1.set_title(title)
                xdata, ydata = choose_window(data, dates, i, windows)
                l.set_ydata(ydata)
                l.set_xdata(xdata)
                l.axes.set_xlim(xdata[0], xdata[len(xdata)-1])
                xdata, ydata = choose_window(data_o, dates, i, windows)
                k.set_ydata(ydata)
                k.set_xdata(xdata)
                k.axes.set_xlim(xdata[0], xdata[len(xdata)-1])
                if autoscale:
                    ax1.set_ylim(ydata.min() - abs(0.1 * ydata.min()),
                                 ydata.max() + 0.1 * ydata.max())
                elif ylim is not None:
                    ax1.set_ylim(ylim[0], ylim[1])
                plt.draw()
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def prev(self, event):
            if len(PickByHand) % 2 == 0:
                self.ind -= 1
                noi = len(windows)-1
                i = self.ind % noi
                if i == 0:
                    msg = 'First timewindow of station: '
                    ax1.set_title(msg + title)
                else:
                    ax1.set_title(title)
                xdata, ydata = choose_window(data, dates, i, windows)
                l.set_ydata(ydata)
                l.set_xdata(xdata)
                l.axes.set_xlim(xdata[0], xdata[len(xdata)-1])
                xdata, ydata = choose_window(data_o, dates, i, windows)
                k.set_ydata(ydata)
                k.set_xdata(xdata)
                k.axes.set_xlim(xdata[0], xdata[len(xdata)-1])
                if autoscale:
                    ax1.set_ylim(ydata.min() - abs(0.1 * ydata.min()),
                                 ydata.max() + 0.1 * ydata.max())
                elif ylim is not None:
                    ax1.set_ylim(ylim[0], ylim[1])
                plt.draw()
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def delete(self, event):
            if len(pick_lines) != 0:
                # if PickByHand[-1] > dates[windows[self.ind]]:
                #     if PickByHand[-1] < dates[windows[self.ind+1]]:
                PickByHand.pop()
                pick_lines[-1].remove()
                pick_lines.pop()
                print('Last pick deleted')
                plt.draw()
            else:
                return

        def quit(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = True
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def nextstation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                tr_i[0] = 1
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def prevstation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                tr_i[0] = -1
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def deletestation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                delete[0] = True
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

    data = trace.data
    data_o = trace_org.data

    fig, ax1 = plt.subplots(1, 1)
    fig.set_size_inches(14, 10)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)

    if ylim:
        ax1.set_ylim(ylim[0], ylim[1])

    if yaxis == 'dates':
        ts_start = trace.stats.starttime.timestamp
        ts_end = trace.stats.endtime.timestamp
        timestamps = np.linspace(ts_start, ts_end, len(data))
        dates = [dt.datetime.utcfromtimestamp(ts) for ts in timestamps]
        ax = plt.gca()
        xfmt = mpld.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.xaxis.set_major_formatter(xfmt)
        plt.xticks(rotation=25)

    else:
        ts_start = trace.stats.starttime
        ts_end = trace.stats.endtime
        t_diff = (ts_end - ts_start) / 3600.
        dates = np.linspace(0, t_diff, len(data))
        ax1.set_xlabel('hours')

    k, = ax1.plot(dates[windows[0]:windows[1]], data_o[windows[0]:windows[1]],
                  color='black', label='%s_org' % trace_org.stats.station)

    l, = ax1.plot(dates[windows[0]:windows[1]], data[windows[0]:windows[1]],
                  picker=5, label='%s_ins' % trace.stats.station)

    ax1.legend()

    global PickByHand, quit_prog
    quit_prog = [False]
    delete = [False]
    tr_i = [1]
    PickByHand = []
    pick_lines = []

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            pick_tmp = zip(np.take(xdata, ind), np.take(ydata, ind))
            amax = np.array(pick_tmp).transpose()[1].argmax()
            pick = pick_tmp[amax]
            ax1.set_ylim(ax1.get_ylim())
            p, = ax1.plot((pick[0], pick[0]), (ax1.get_ylim()), color='red')
            pick_lines.append(p)
            plt.draw()

            # Add confirmation here "Keep pick?"
            # Add feedback in plot (Mark lines of pick)
            # count 2 clicks at then ask
            PickByHand.append(pick[0])

    fig.canvas.mpl_connect('pick_event', onpick1)
    callback = Index()

    axquit = plt.axes([0.91, 0.9, 0.08, 0.075])

    axnext = plt.axes([0.91, 0.75, 0.08, 0.075])
    axprev = plt.axes([0.91, 0.65, 0.08, 0.075])

    axdelete = plt.axes([0.91, 0.55, 0.08, 0.075])

    axnextstat = plt.axes([0.91, 0.35, 0.08, 0.075])
    axprevstat = plt.axes([0.91, 0.25, 0.08, 0.075])
    axremstat = plt.axes([0.91, 0.15, 0.08, 0.075])

    bnext = Button(axnext, 'Next TW')
    bnext.on_clicked(callback.next)

    bprev = Button(axprev, 'Previous TW')
    bprev.on_clicked(callback.prev)

    bdelete = Button(axdelete, 'Remove Pick')
    bdelete.on_clicked(callback.delete)

    bnextstat = Button(axnextstat, 'Next Station')
    bnextstat.on_clicked(callback.nextstation)

    bprevstat = Button(axprevstat, 'Previous Station')
    bprevstat.on_clicked(callback.prevstation)

    bquit = Button(axquit, 'Quit')
    bquit.on_clicked(callback.quit)

    brem = Button(axremstat, 'Delete Station')
    brem.on_clicked(callback.deletestation)

    plt.show()

    return PickByHand, tr_i[0], quit_prog[0], delete[0]


def pick_inspector_window(traces, windows, xlabel, ylabel, title,
                          ylim=None, autoscale=False, cursor=False):
    global PickByHand, quit_prog
    fstart = 100
    quit_prog = [False]
    delete = [False]
    tr_i = [1]
    PickByHand = []
    pick_lines = []

    msg = np.nditer([['Last pick deleted -',
                      'Last pick deleted \\',
                      'Last pick deleted |',
                      'Last pick deleted /']])

    def trcopy(x):
        y = x.copy()
        return y

    def choose_window(data, dates, i, windows):
        win_trace = data[windows[i]:windows[i+1]]
        win_dates = dates[windows[i]:windows[i+1]]
        return win_dates, win_trace

    def _mask_data(ind, tr_work):
        for win in pairwise(PickByHand):
            for i, t in enumerate(tr_work):
                i1 = int(np.round(win[0]*3600. * t.stats.sampling_rate))
                i2 = int(np.round(win[1]*3600. * t.stats.sampling_rate))
                t.data = mask_data(t.data, i1, i2, shape='linear')

        for l, t in zip(lines, traces):
            xdata, ydata = choose_window(t.data, dates, ind, windows)
            l.set_ydata(ydata)
            l.set_xdata(xdata)

        for l, t in zip(ftlines, traces):
            ft = fourier_transform(t, tw[0], tw[1], 0, nend, 'hanning')[1]
            l.set_ydata(abs(ft[fstart:len(ft)/2]))

    class Index(object):
        ind = 0
        trbak = trcopy(traces)

        # Side panel ##############################################
        def quit(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = True
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def nextstation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                tr_i[0] = 1
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def prevstation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                tr_i[0] = -1
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def deletestation(self, event):
            if len(PickByHand) % 2 == 0:
                quit_prog[0] = False
                delete[0] = True
                plt.close('all')
            else:
                print('Uneven window, picks must be a mutliple of 2')

        def undo(self, event):
            i = self.ind
            data = []
            for t, d in zip(traces, self.trbak):
                t.data = d.data.copy()
                data.append(d.data)

            for l, d in zip(lines, data):
                xdata, ydata = choose_window(d, dates, i, windows)
                l.set_ydata(ydata)
                l.set_xdata(xdata)

            for l, t in zip(ftlines, traces):
                ft = fourier_transform(t, tw[0], tw[1], 0, nend, 'hanning')[1]
                l.set_ydata(abs(ft[fstart:len(ft)/2]))

        # Bottom Panel ##############################################
        def ecomp(self, event):
            self.trbak = traces.copy()
            tr_work = traces.select(channel='*E')
            tr_work = _mask_data(self.ind, tr_work)

        def ncomp(self, event):
            self.trbak = traces.copy()
            tr_work = traces.select(channel='*N')
            tr_work = _mask_data(self.ind, tr_work)

        def zcomp(self, event):
            self.trbak = traces.copy()
            tr_work = traces.select(channel='*Z')
            tr_work = _mask_data(self.ind, tr_work)

        def hcomp(self, event):
            self.trbak = traces.copy()
            tr_work = traces.select(channel='*E')
            tr_work += traces.select(channel='*N')
            tr_work = _mask_data(self.ind, tr_work)

        def allcomp(self, event):
            self.trbak = traces.copy()
            _mask_data(self.ind, traces)

        def deletepicks(self, event):
            if len(pick_lines) != 0:
                # if PickByHand[-1] > dates[windows[self.ind]]:
                #     if PickByHand[-1] < dates[windows[self.ind+1]]:
                PickByHand.pop()
                pick_lines[-1].remove()
                pick_lines.pop()
                try:
                    print('%s' % next(msg), end='\r')
                    sys.stdout.flush()
                except StopIteration:
                    msg.reset()
                    print('%s' % next(msg), end='\r')
                    sys.stdout.flush()
                plt.draw()

        def detrend(self, event):
            self.trbak = traces.copy()
            i = self.ind
            traces.detrend('constant')
            traces.detrend('linear')
            data = []
            for t in traces:
                data.append(t.data)
            for l, d in zip(lines, data):
                xdata, ydata = choose_window(d, dates, i, windows)
                l.set_ydata(ydata)
                l.set_xdata(xdata)

    callback = Index()

    ts_start = traces[0].stats.starttime
    ts_end = traces[0].stats.endtime
    t_diff = (ts_end - ts_start) / 3600.
    dates = np.linspace(0, t_diff, len(traces[0].data))

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2)
    ax = range(2)
    ax[0] = plt.subplot(gs[0, 0])
    ax[1] = plt.subplot(gs[0, 1])

    plt.subplots_adjust(bottom=0.2)
    fig.set_size_inches(14, 10)
    ax[0].set_title(title)
    ax[0].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel)

    ax[1].set_xlabel('Frequency (mHz)')
    if ylim:
        ax[0].set_ylim(ylim[0], ylim[1])

    ax[0].set_xlabel('hours')

    lines = []
    ftlines = []
    data = []
    chlabel = []
    FFT = []
    f = []
    for t in traces:
        data.append(t.data)
        tw = [0, (t.stats.endtime-t.stats.starttime)/3600.]
        nend = t.stats.npts
        freq, Fxx, delomeg = fourier_transform(t, tw[0], tw[1],
                                               0, nend, 'hanning')
        FFT.append(Fxx[fstart:len(Fxx)/2])
        f.append(freq[:len(f)/2])
        label = '%s' % (t.stats.channel)

        tsl = ax[0].plot(dates[windows[0]:windows[1]],
                         t.data[windows[0]:windows[1]],
                         picker=1, label=label)
        lines.append(tsl[0])
        ax[0].set_xlim(t.data[0], t.data[len(t.data)-1])
        chlabel.append(label)
        xdata, ydata = choose_window(t.data, dates, 0, windows)
        ax[0].set_xlim(xdata[0], xdata[len(xdata)-1])

        ftl = ax[1].plot(freq[fstart:len(Fxx)/2],
                         abs(Fxx[fstart:len(Fxx)/2]), label=label)
        ftlines.append(ftl[0])

    ax[0].legend()
    ax[1].legend()
    checkbuttons = plt.axes([0.91, 0.5, 0.08, 0.075])
    check = CheckButtons(checkbuttons, chlabel, (True, True, True))

    def func(label):
        for ch, l in zip(chlabel, lines):
            if label == ch:
                l.set_visible(not l.get_visible())
        for ch, l in zip(chlabel, ftlines):
            if label == ch:
                l.set_visible(not l.get_visible())
        plt.draw()

    check.on_clicked(func)

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            if len(ind) > 1:
                datax, datay = event.artist.get_data()
                datax, datay = [datax[i] for i in ind], [datay[i] for i in ind]
                msx, msy = event.mouseevent.xdata, event.mouseevent.ydata
                dist = np.sqrt((np.array(datax)-msx)**2
                               + (np.array(datay)-msy)**2)
                ind = [ind[np.argmin(dist)]]
            pick_tmp = zip(np.take(xdata, ind), np.take(ydata, ind))
            amax = np.array(pick_tmp).transpose()[1].argmax()
            pick = pick_tmp[amax]
            ax[0].set_ylim(ax[0].get_ylim())
            p, = ax[0].plot((pick[0], pick[0]), (ax[0].get_ylim()),
                            color='red')
            pick_lines.append(p)
            plt.draw()

            # Add confirmation here "Keep pick?"
            # Add feedback in plot (Mark lines of pick)
            # count 2 clicks at then ask
            PickByHand.append(pick[0])
            PickByHand.sort()

    fig.canvas.mpl_connect('pick_event', onpick1)

    # Side Panel ##############################################
    axquit = plt.axes([0.11, 0.9, 0.08, 0.05])
    # axtwnext = plt.axes([0.91, 0.75, 0.08, 0.075])
    # axtwprev = plt.axes([0.91, 0.65, 0.08, 0.075])
    axnextstat = plt.axes([0.91, 0.35, 0.08, 0.075])
    axprevstat = plt.axes([0.91, 0.25, 0.08, 0.075])
    axundo = plt.axes([0.91, 0.15, 0.08, 0.075])
    axdeletestat = plt.axes([0.81, 0.9, 0.08, 0.075])

    # Bottom Panel ##############################################

    epick = plt.axes([0.11, 0.01, 0.08, 0.05])
    npick = plt.axes([0.21, 0.01, 0.08, 0.05])
    zpick = plt.axes([0.31, 0.01, 0.08, 0.05])
    horpick = plt.axes([0.41, 0.01, 0.08, 0.05])
    allpick = plt.axes([0.51, 0.01, 0.08, 0.05])
    axdetrend = plt.axes([0.61, 0.01, 0.08, 0.05])
    rmpick = plt.axes([0.81, 0.01, 0.08, 0.05])

    # Callbacks ###############################################
    bquit = Button(axquit, 'Quit')
    bquit.on_clicked(callback.quit)

    # bnext = Button(axtwnext, 'Next TW')
    # bnext.on_clicked(callback.next)
    #
    # bprev = Button(axtwprev, 'Previous TW')
    # bprev.on_clicked(callback.prev)

    bnextstat = Button(axnextstat, 'Next Station')
    bnextstat.on_clicked(callback.nextstation)

    bprevstat = Button(axprevstat, 'Prev. Station')
    bprevstat.on_clicked(callback.prevstation)

    brem = Button(axdeletestat, 'Delete Station')
    brem.on_clicked(callback.deletestation)

    bundo = Button(axundo, 'Undo')
    bundo.on_clicked(callback.undo)

    bepick = Button(epick, 'E')
    bepick.on_clicked(callback.ecomp)

    bnpick = Button(npick, 'N')
    bnpick.on_clicked(callback.ncomp)

    bzpick = Button(zpick, 'Z')
    bzpick.on_clicked(callback.zcomp)

    bhpick = Button(horpick, 'Hor.')
    bhpick.on_clicked(callback.hcomp)

    ballpick = Button(allpick, 'All')
    ballpick.on_clicked(callback.allcomp)

    bdetrend = Button(axdetrend, 'Detrend')
    bdetrend.on_clicked(callback.detrend)

    brmpick = Button(rmpick, 'Remove Pick')
    brmpick.on_clicked(callback.deletepicks)

    plt.show()
    return traces, PickByHand, tr_i[0], quit_prog[0], delete[0]
