try:
    import cartopy
    import cartopy.crs as ccrs
    from frospy.plot.maps import plot_inv_cat_c2py, gcpmap_c2py, gcp2map
except Exception as e:
    print("Error in cartopy, not imported!\n%s" % e)

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker

from frospy.plot.nmplt import (axvlines, axtexts, format_exponent,
                             get_text_positions, text_plotter)
from frospy.util.read import read_std_cat

import numpy as np
from obspy.core import AttribDict

from frospy.core import Spectrum


def plot_spectrum(main, gui=False):
    """
    Plot routine for spectrum
    cmap: options: 'rainbow' or 'random'
    """
    fig = main.rfig
    ax_old = main.rax
    seg_ax = main.seg_ax

    if not fig and gui is False:
        plt.ion()

    # Setting title
    title = main.spec.stats.__str__()
    title.encode('utf8')

    # Getting Parameters for current plot
    slat = round(main.spec.stats.record.ah.station.latitude, 2)
    slon = round(main.spec.stats.record.ah.station.longitude, 2)
    elat = round(main.spec.stats.record.ah.event.latitude, 2)
    elon = round(main.spec.stats.record.ah.event.longitude, 2)

    prox = 15.
    if main.minispec == 'amp':
        gs = None
        n_plots = 1
    elif main.minispec is True:
        gs = gridspec.GridSpec(6, 1)
        n_plots = 2
    elif main.minispec == 'lat':
        gs = gridspec.GridSpec(12, 1)
        n_plots = 2
    elif main.minispec == 'pretty':
        gs = gridspec.GridSpec(3, 1)
        n_plots = 2
    else:
        gs = gridspec.GridSpec(12, 6)
        n_plots = 9

    if fig and plt.fignum_exists(fig.number):
        plt.figure(fig.number)
        for no, axis in enumerate(fig.axes[:n_plots]):
            axis.clear()
        ax = ax_old

    else:
        fig, ax, seg_ax = init_plot_spectrum(gs, main.minispec)
        if main.minispec != 'pretty':
            fig.set_size_inches(main.fig_size)

    if main.minispec is False:
        ax.real.set_title(title)
    elif main.minispec != 'pretty':
        title += "\n (%s / %s)\n" % (main.i+1, len(main.st_work))
        ax.amp.set_title(title)

    if main.minispec is True:
        # Plotting Phase
        main.spec.plot(main.fw[0], main.fw[1], part='Phase',
                       width=main.line_width, ax=ax.phase,
                       cmap=main.cmap, cmap_highlight=main.cmap_highlight)
    if main.minispec == 'pretty':
        main.spec.plot(main.fw[0], main.fw[1], part='Phase',
                       width=main.line_width, ax=ax.phase,
                       cmap=main.cmap, xlabel='f(mHz)', normalize=True,
                       cmap_highlight=main.cmap_highlight)

    if main.minispec is False:
        # Plotting Phase
        main.spec.plot(main.fw[0], main.fw[1], part='Phase',
                       width=main.line_width, ax=ax.phase,
                       cmap=main.cmap, cmap_highlight=main.cmap_highlight)
        # Plotting Real part
        main.spec.plot(main.fw[0], main.fw[1], part='Re', width=main.line_width,
                       ax=ax.real, cmap=main.cmap,
                       cmap_highlight=main.cmap_highlight)

        # Plotting Imaginary part
        main.spec.plot(main.fw[0], main.fw[1], part='Im', width=main.line_width,
                       ax=ax.imag, cmap=main.cmap,
                       cmap_highlight=main.cmap_highlight)

    # Plotting Amplitude
    dlabel = spectrum_label(main.spec)
    if main.minispec == 'pretty':
        main.spec.plot(main.fw[0], main.fw[1], part='Amplitude',
                       xlabel='f(mHz)',
                       ylabel='Amplitude', normalize=True,
                       ax=ax.amp, width=main.line_width,
                       cmap=main.cmap, cmap_highlight=main.cmap_highlight)
    else:
        main.spec.plot(main.fw[0], main.fw[1], part='Amplitude',
                       xlabel='frequency (mHz)',
                       ylabel='Normalized Amplitude',
                       ax=ax.amp, width=main.line_width, dlabel=dlabel,
                       cmap=main.cmap,
                       cmap_highlight=main.cmap_highlight)

    # Plotting Amplitude: Lines of modes
    if main.modes:
        if main.minispec is not False:
            l_height = 0.11
            l_width = 0.012
        else:
            l_height = 0.06
            l_width = 0.03
        if main.minispec == 'pretty':
            ypos = 0.4
            for mode in main.modes:
                mtxt = '${}_{%d}%s_{%d}$' % (mode.n, mode.type, mode.l)
                fig.axes[1].text(mode.freq, ypos, mtxt, fontsize=14)
                ypos += 0.07
        else:
            plot_modes(main.spec, main.fw[0], main.fw[1], main.modes, ax.amp,
                       l_height, l_width, fontsize=main.fs)

    if main.segments:
        plot_segments(main.spec, main.segments, main.fw[0], main.fw[1], seg_ax)

    if main.noisewin:
        plot_noise_window(main.spec, main.fw[0], main.fw[1], ax.amp)

    ax.amp.legend()

    if main.minispec is not True:
        if main.minispec not in ('amp', 'lat', 'pretty'):
            # Plotting seismogram
            plot_seismogram(main.spec, ax.seis)

    for axis in seg_ax:
        axis = format_exponent(axis)

    # Plotting map with GCP and CMT
    station = main.spec.stats.station.code

    # Plot Station map
    if main.minispec is False:
        if main.inv is not None:
            ax['stat'] = init_map(gs[0:4, 2:4])
            plot_station_map(main.spec, main.inv, main.segments, station,
                             gs, ax.stat)
        # else:
        #     ax.stat.axis('off')

        ax['gcp'] = init_gcpmap(gs[4:8, 2:4], slon, elon)
        plot_gcp_map(slat, slon, elat, elon, main.cmt, main.cmt_id,
                     prox, ax.gcp)

        # Magnitude information
        plot_magnitude(main.cmt, ax.mag)

    if main.minispec is False or main.minispec == 'lat':
        # Latitude Spectrum plot
        plot_speclat(main.st, main.spec, [main.fw[0], main.fw[1]],
                     zoom=main.zoom_sl, segments=main.segments, highlight=None,
                     modes=main.modes, ax=ax.speclat, print_labels=False,
                     gui=gui)

    # fig.tight_layout()
    if main.minispec == 'pretty':
        fig.subplots_adjust(hspace=0)
        fig.axes[0].tick_params(axis='both', which='both',
                                bottom=False, labelsize=main.fs - 2, width=1.4)
        fig.axes[0].spines['bottom'].set_visible(False)
        for axis in ['top', 'bottom', 'left', 'right']:
            fig.axes[0].spines[axis].set_linewidth(1.8)
            fig.axes[1].spines[axis].set_linewidth(1.8)
        fig.axes[1].get_legend().remove()
        fig.axes[1].xaxis.set_major_locator(plt.MaxNLocator(3))
        fig.axes[1].tick_params(axis='both', which='major',
                                labelsize=main.fs - 2, width=1.4)
        fig.axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
        fig.axes[1].set_xlabel("f(mHz)", fontsize=main.fs)
        fig.axes[0].set_ylabel("Phase", fontsize=main.fs)
        fig.axes[1].set_ylabel("Amplitude", fontsize=main.fs)
        cat = read_std_cat(main.cmt_id)

        for descr in cat[0].event_descriptions:
            if descr['type'] == 'region name':
                evloc = descr['text']
        yr = main.cmt_id[4:6]
        if yr.startswith(('7', '8', '9')):
            yr = '19'+yr
        else:
            yr = '20'+yr
        t1 = ''
        if main.modes:
            if len(main.modes) > 1:
                t1 = 'Modes'
            else:
                t1 = 'Mode'
            for i, mode in enumerate(main.modes):
                if i == 1:
                    t1 += " ${}_{%d}%s_{%d}$" % (mode.n, mode.type, mode.l)
                else:
                    t1 += "-${}_{%d}%s_{%d}$" % (mode.n, mode.type, mode.l)
        t2 = "$M_w$ $%.1f,$ %s$, %s$" % (main.cmt[6], evloc, yr)
        t3 = "$%s$ $\mathrm{station,}$ $%d$-$%d\mathrm{h}$"
        t3 = t3 % (main.spec.stats.station.code, main.spec.stats.tw[0]/3600.,
                   main.spec.stats.tw[1]/3600.)
        if main.fig_abc is not False:
            pretty_title = r"$\bf{" + main.fig_abc + ")}$  %s\n%s" % (t2, t3)
        else:
            pretty_title = '%s\n%s\n%s' % (t1, t2, t3)

        fig.axes[0].set_title(pretty_title, fontsize=main.fs)

        fig.axes[1].legend(loc='upper left', frameon=False, borderpad=0,
                           handlelength=1, handletextpad=0.3, labelspacing=0.2,
                           fontsize=main.fs - 2)
        fig.set_size_inches(5.5, 6.5)
        plt.subplots_adjust(left=0.1, top=0.82)

    if main.savefig:
        msg = "\033[92mSaving plotfile..\033[0m"
        print(msg)
        fname = '%s-%s.%s.png' % (main.cmt_id, main.spec.stats.station.code,
                                  main.spec.stats.record.channel)
        plt.ioff()
        # fig.set_size_inches(20, 12)
        fig.savefig(fname, dpi=300, orientation='portrait')
        plt.close(fig)
        fig = None
        ax = None
    else:
        fig.canvas.draw()
        if gui is False:
            plt.show()

    mpl.rcParams.update({'font.size': main.fs})
    for _ax in ax.values():
        _ax.tick_params(axis='both', which='major', labelsize=main.fs)
        [i.set_linewidth(main.border_width) for i in _ax.spines.values()]
        _ax.xaxis.set_tick_params(which='both', width=main.tick_width)
        _ax.yaxis.set_tick_params(which='both', width=main.tick_width)
        _ax.xaxis.label.set_size(main.fs)
        _ax.yaxis.label.set_size(main.fs)

    if main.ylim is not None:
        ax.amp.set_ylim(main.ylim)
    main.rfig = fig
    main.rax = ax
    main.seg_ax = seg_ax

    return main


def init_plot_spectrum(gs, minispec):
    ax = AttribDict()

    if minispec == 'amp':
        fig, ax['amp'] = plt.subplots()
        seg_ax = [ax.amp]

    elif minispec == 'lat':
        fig = plt.figure()
        ax['speclat'] = plt.subplot(gs[6:12, 0])
        ax['amp'] = plt.subplot(gs[0:5, 0])
        seg_ax = [ax.amp]

    elif minispec is True:
        fig = plt.figure()
        # Phase spectrum
        ax['phase'] = plt.subplot(gs[0:2, 0])
        # # Timeseries
        # ax['seis'] = plt.subplot(gs[0:2, 0])
        # Amplitudes
        ax['amp'] = plt.subplot(gs[3:5, 0])
        seg_ax = [ax.amp, ax.phase]

    elif minispec == 'pretty':
        fig, (ax['phase'], ax['amp']) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios':[1,3]})
        seg_ax = [ax.amp, ax.phase]


    elif minispec == 'gui':
        fig, (ax['phase'], ax['amp']) = plt.subplots(2, sharex=True)
        seg_ax = [ax.amp, ax.phase]

    else:
        fig = plt.figure()
        # # Real part
        ax['real'] = plt.subplot(gs[0:2, 0:2])

        # # Imaginary part
        ax['imag'] = plt.subplot(gs[2:4, 0:2])

        # Phase
        ax['phase'] = plt.subplot(gs[4:6, 0:2])

        # Amplitude
        ax['amp'] = plt.subplot(gs[6:12, 0:2])

        # Timeseries
        ax['seis'] = plt.subplot(gs[10:12, 2:4])

        # Magnitudes
        ax['mag'] = plt.subplot(gs[8:10, 2:4])

        # These axis is set in maps
        # ax[4] = init_gcpmap(gs[4:8, 2:4], slon, elon)
        #
        # if inv:
        #     ax[5] = init_map(gs[0:4, 2:4])
        # Latitude plot of main.st_work
        ax['speclat'] = plt.subplot(gs[0:12, 4:6])
        seg_ax = [ax.real, ax.imag, ax.phase, ax.amp]
    return fig, ax, seg_ax


def init_map(gs):
    """
    param grd_size: defining the grid size
    type  grd_size: 2x2 tupel

    param grd_pos: position on grid
    type  grd_pos: 2x2 tupel
    """
    ocean_color = '#EBEBEB'
    land_color = '#FBFBF2'
    if gs is not None:
        ax = plt.subplot(gs, projection=ccrs.Mollweide())
    else:
        ax = plt.subplot(projection=ccrs.Mollweide())
    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()
    # ax.stock_img()
    ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    ax.coastlines()
    ax.gridlines()
    return ax


def init_gcpmap(gs, lon1, lon2):
    ocean_color = '#EBEBEB'
    land_color = '#FBFBF2'

    difflon = abs(lon1 - lon2)
    if difflon > 180:
        lon0 = int((lon1+lon2)/2. + 180)
    else:
        lon0 = int((lon1+lon2)/2.)
    proj = ccrs.Mollweide(central_longitude=lon0)

    if gs is not None:
        ax = plt.subplot(gs, projection=proj)
    else:
        ax = plt.subplot(projection=proj)
    ax.set_global()
    # ax.stock_img()
    ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    ax.coastlines()
    ax.gridlines()
    return ax


def plot_modes(spectrum, fw1, fw2, modes, ax, l_height, l_width, fontsize=10):

    startlabel = spectrum.flabel(fw1)
    endlabel = spectrum.flabel(fw2)
    f = spectrum.stats.freq
    mode_lines(ax, f[startlabel:endlabel+1], modes, label_height=l_height,
               label_width=l_width, fontsize=fontsize)


def plot_segments(spectrum, segments, fw1, fw2, ax, alpha=0.05,
                  spec_reference=True):
    startlabel = spectrum.flabel(fw1)
    endlabel = spectrum.flabel(fw2)
    f = spectrum.stats.freq

    if spec_reference is True:
        segments = segments.select(station=spectrum.stats.station.code)
    for pick in segments:
        if pick.fw1 > f[startlabel] and pick.fw2 < f[endlabel+1]:
            for subax in ax:
                subax.axvline(x=pick.fw1, linestyle='-',
                              linewidth=1, color='blue')
                subax.axvline(x=pick.fw2, linestyle='-',
                              linewidth=1, color='blue')
                subax.axvspan(pick.fw1, pick.fw2, alpha=alpha, color='blue')


def plot_noise_window(spectrum, segments, fw1, fw2, ax):
    startlabel = spectrum.flabel(fw1)
    endlabel = spectrum.flabel(fw2)
    f = spectrum.stats.freq
    for pick in segments.select(station=spectrum.stats.station.code):
        nwin = spectrum.signal2noise(pick)[1]
        for i, n in enumerate(nwin):
            if n[0] > f[startlabel] and n[1] < f[endlabel+1]:
                if i == 0:
                    ax.axvspan(n[0], n[1], alpha=0.5, color='red',
                               label='noise')
                else:
                    ax.axvspan(n[0], n[1], alpha=0.5, color='red')


def plot_seismogram(spectrum, ax):
    if spectrum.syn:
        for xsyn in spectrum.syn:
            y = list(xsyn.timeseries.values())[0]
            ax.plot(spectrum.stats.times / 3600., y, linestyle='dashed')

    y = list(spectrum.data.timeseries.values())[0]
    ax.plot(spectrum.stats.times / 3600., y, color='black')
    ax.set_xlabel('time (h)')
    ax.axes.yaxis.set_ticklabels([])

    t_total = (spectrum.stats.record.endtime -
               spectrum.stats.record.starttime) / 3600
    ax.set_title('%i hours of data in Seismogram' % t_total)


def plot_speclat(stream, spectrum, fw, zoom=1, segments=None,
                 highlight=None, modes=None, ax=None, print_labels=True,
                 gui=False,
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
    tw = spectrum.stats.tw / 3600.
    if zoom == 'normalize each':
        max_amp = 0.0
        for tr in stream:
            if segments:
                seg = segments.select(station=tr.stats.station)[0]
                tw = [seg.tw1, seg.tw2]

                if fw is None:
                    fw = [seg.fw1, seg.fw2]
            spec_tmp = Spectrum(tr, tw[0], tw[1])
            amp = abs(spec_tmp.data.fft.Data)
            startlabel = spec_tmp.flabel(fw[0])
            endlabel = spec_tmp.flabel(fw[1])
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

        spec = Spectrum(tr, tw[0], tw[1])
        Fxx = spec.data.fft.Data
        amp = abs(Fxx)

        lat = tr.stats.ah.station.latitude
        f = spectrum.stats.freq
        startlabel = spectrum.flabel(fw[0])
        endlabel = spectrum.flabel(fw[1])

        if i == 0:
            x_extend = [startlabel, endlabel]
        else:
            if startlabel < x_extend[0]:
                x_extend[0] = startlabel
            if endlabel > x_extend[1]:
                x_extend[1] = endlabel

        freq = np.append(f[startlabel:endlabel+1], np.nan)
        amp_tr = amp[startlabel:endlabel+1]
        if type(zoom) == float or type(float) == int:
            amp = np.append(amp_tr*float(zoom)/amp_tr.max() + lat, np.nan)
        else:
            amp = np.append(zoom*amp_tr + lat, np.nan)

        if (
           tr.stats.station in highlight or
           tr.stats.station == spectrum.stats.station.code
           ):
            if freqs_hl is None:
                freqs_hl = freq
                amps_hl = amp
            # else:
            #     freqs_hl = np.vstack((freqs, freq))
            #     amps_hl = np.vstack((amps, amp))
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
    if modes is not None:
        mode_lines(ax, f[startlabel:endlabel+1], modes, label_height=0.025,
                   label_width=0.03, print_labels=print_labels)
    if segments:
        plot_segments(spectrum, segments, fw[0], fw[1], [ax])

    # Loop for station names
    ytrans = ax.get_yaxis_transform()
    for station, lat in lat_label.items():
        ax.annotate("%s" % station, xy=(0.05, lat), xycoords=ytrans)

    ax.set_xlim(f[x_extend[0]], f[x_extend[1]+1])
    if gui is False:
        plt.show()
    return


def spectrum_label(spectrum):
    # If default label is set, in spectrum, print default label in plot
    if 'data' not in spectrum.data.fft.keys():
        d = list(spectrum.data.fft.keys())[0]
    else:
        d = '%s.%s.%s.%s' % (spectrum.stats.record.network,
                             spectrum.stats.record.station,
                             spectrum.stats.record.location,
                             spectrum.stats.record.channel)
    return d


def mode_lines(ax, xrange, modes, overlap=True,
               label_height=None, label_width=None, print_labels=True,
               fontsize=10):
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
    if mode_labels and print_labels is True:
        mode_labels = np.array(mode_labels, dtype=object)
        m_names = mode_labels.transpose()[1]

        if label_height is None or label_width is None:
            axtexts(mode_freqs, m_names,
                    y=1.03, ax=ax, transform=xtrans, rotation=45,
                    ha="center", va='center')
        else:
            m_ypos = (np.ones(len(mode_freqs)) * ax.get_ylim()[1])
            txt_height = label_height * (ax.get_ylim()[1] - ax.get_ylim()[0])
            txt_width = 1.2 * label_width * (ax.get_xlim()[1] - ax.get_xlim()[0])
            text_positions = get_text_positions(mode_freqs, m_ypos, txt_width,
                                                txt_height)
            text_plotter(mode_freqs, m_ypos, text_positions, m_names, ax,
                         txt_width, txt_height, fontsize)

    return


def plot_station_map(spectrum, inv, segments, station, gs, ax):
    if segments:
        no_used_sta = len(segments)
        ax = plot_inv_cat_c2py(inv, None, ax, station, segments)
        ax.set_title('%s/%s Stations' % (no_used_sta, len(inv[0])))
        ax.axis('off')
    else:
        ax = plot_inv_cat_c2py(inv, None, ax, station, None)
        ax.set_title('%s Stations' % len(inv[0]))
        ax.axis('off')


def plot_gcp_map(slat, slon, elat, elon, cmt, cmt_id, prox, ax,
                 basemap=False):
    if basemap:
        ax = gcp2map(slat, slon, elat, elon, ax, cmt=cmt,
                     projection='moll', proximity=prox)
    else:
        ax = gcpmap_c2py(slat, slon, elat, elon, ax, cmt=cmt,
                         proximity=prox)
    ax.set_title('GCP with Event: %s' % cmt_id)
    ax.axis('off')


def plot_magnitude(cmt, ax):
    ax.axis('off')
    if cmt and len(cmt) >= 6:
        ax.set_title('Centroid Moment Tensor')
        ax.annotate(r"$M_W$: %.1f" % cmt[6], xy=(0.1, 0.8))
        ax.annotate(r"$M_0$: %.2E dyne-cm" % cmt[7], xy=(0.1, 0.4))
        ax.annotate(r"Depth: %.2f km" % cmt[8], xy=(0.1, 0.0))


def plot_spectrum_partials(main):
    """
    Plot routine for spectrum
    cmap: options: 'rainbow' or 'random'
    """

    # Temporary Solution with arguments
    spectrum = main.spec
    fw = main.fw
    modes = main.modes
    segments = main.segments
    st_work = main.st_work
    inv = main.inv
    cmt = main.cmt
    cmt_id = main.cmt_id
    cmap = main.cmap
    fig = main.fig
    seg_ax = main.seg_ax
    minispec = main.minispec

    mpl.rcParams.update({'font.size': main.fs})
    if not fig:
        plt.ion()

    # Setting title
    title = spectrum.stats.__str__()
    title.encode('utf8')

    # Getting Parameters for current plot
    slat = round(spectrum.stats.record.ah.station.latitude, 2)
    slon = round(spectrum.stats.record.ah.station.longitude, 2)
    elat = round(spectrum.stats.record.ah.event.latitude, 2)
    elon = round(spectrum.stats.record.ah.event.longitude, 2)

    prox = 15.
    l_height = 0.06
    l_width = 0.03
    seg_ax = ['real', 'imag', 'phase', 'amp']

    for part in minispec:
        if part == 'stat':
            station = spectrum.stats.station.code
            fig = plt.figure()
            ax = init_map(None)
            plot_station_map(spectrum, inv, segments, station, None, ax)

        elif part == 'gcp':
            fig, ax = plt.subplots()
            basemap = True
            plot_gcp_map(slat, slon, elat, elon, cmt[0:6], cmt_id, prox, ax,
                         basemap)

        else:
            fig, ax = plt.subplots()
            if part == 'amp':
                dlabel = spectrum_label(spectrum)
                spectrum.plot(fw[0], fw[1], part='Amplitude',
                              xlabel='frequency (mHz)',
                              ylabel='Amplitude', normalize=True,
                              ax=ax, width=main.line_width, dlabel=dlabel,
                              cmap=cmap)
                if modes is not None:
                    plot_modes(spectrum, fw[0], fw[1], modes, ax,
                               l_height, l_width, fontsize=main.fs)
                if segments:
                    plot_segments(spectrum, segments, fw[0], fw[1], [ax])
                ax.legend()

            elif part == 'phase':
                dlabel = spectrum_label(spectrum)
                spectrum.plot(fw[0], fw[1], part='Phase', width=main.line_width,
                              ax=ax, cmap=cmap, normalize=True)
                if modes is not None:
                    plot_modes(spectrum, fw[0], fw[1], modes, ax,
                               l_height, l_width, fontsize=main.fs)
                if segments:
                    plot_segments(spectrum, segments, fw[0], fw[1], [ax])

            elif part == 'mag':
                print(cmt)
                plot_magnitude(cmt, ax)

            elif part == 'lat':
                plot_speclat(st_work, spectrum, [fw[0], fw[1]],
                             zoom=main.zoom_sl,
                             segments=segments, highlight=None,
                             modes=modes, ax=ax)

            else:
                spectrum.plot(fw[0], fw[1], part=part, width=main.line_width,
                              ax=ax, cmap=cmap)

            if part in seg_ax:
                ax = format_exponent(ax)
        fig.tight_layout()

    plt.show()

    main.rfig = fig
    main.rax = ax
    main.seg_ax = seg_ax
    return main
