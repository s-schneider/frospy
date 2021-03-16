from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
import matplotlib.pyplot as plt

try:
    import cartopy
    import cartopy.crs as ccrs
    from cartopy.io.img_tiles import GoogleTiles
except ImportError:
    CTP = False

from matplotlib.pyplot import cm

from frospy.util.array_util import geometrical_center
from frospy.plot.nmplt import beachball
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
import numpy as np

# colors
global ocean_color
global land_color
global all_station_c
global used_station_c
global unused_station_c
global marked_station_c
global event_color
global marker_size

ocean_color = '#EBEBEB'
land_color = '#FBFBF2'
# all_station_c = '#2A98FF'
used_station_c = '#e57248'
unused_station_c = '#bc8f8f'
all_station_c = '#e57248'
# used_station_c = '#e57248'
# unused_station_c = '#bc8f8f'
marked_station_c = 'yellow'
event_color = '#'
marker_size = 4


class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/'
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url


def plot_gcp_old(inv, event, phaselist=None, level='array'):
    # lon_0 is central longitude of projection, lat_0 the central latitude.
    # resolution = 'c' means use crude resolution coastlines, 'l' means low,
    # 'h' high etc.
    # zorder is the plotting level, 0 is the lowest, 1 = one level higher ...
    # m = Basemap(projection='nsper',lon_0=20, lat_0=25,resolution='c')
    m = Basemap(projection='kav7', lon_0=0, resolution='c')
    m.drawmapboundary(fill_color='#B4FFFF')
    m.fillcontinents(color='#00CC00', lake_color='#B4FFFF', zorder=0)

    tm = TauPyModel('ak135')
    elat = event.origins[0].latitude
    elon = event.origins[0].longitude
    depth = event.origins[0].depth / 1000.
    for net in inv:
        if level in ['array']:
            geoc = geometrical_center(net)
            rlat = [geoc.latitude]
            rlon = [geoc.latitude]

            if phaselist:
                pp = tm.get_pierce_points_geo(depth, elat, elon, rlat, rlon,
                                              phaselist)
                px = pp.pierce
                py = pp.pierce
                m.scatter(px, py, 10, marker='d', color='yellow', zorder=2)
        elif level in ['station']:
            rlat = []
            rlon = []
            for station in net:
                rlat.append(station.latitude)
                rlon.append(station.longitude)

        for lats, lons in zip(rlat, rlon):
            # import station coordinates, with symbol (^ = triangle)
            m.scatter(lats, lons, 800, marker='^', color='#004BCB', zorder=1)
            m.drawgreatcircle(elon, elat, lons, lats, linewidth=1,
                              color='black', zorder=1)

        # import event coordinates, with symbol (* = Star)
        m.scatter(elat, elon, 800, marker='*', color='red', zorder=1)

    m.drawcoastlines(zorder=1)

    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 30.), zorder=1)
    m.drawmeridians(np.arange(0., 420., 60.), zorder=1)
    plt.title("")
    plt.show()


def plot_gcp(inv, event, phaselist=None, level='station', map=None,
             proj='kav7'):

    elat = event.origins[0].latitude
    elon = event.origins[0].longitude
    for net in inv:
        if level in ['array']:
            geoc = geometrical_center(net)
            slat = [geoc.latitude]
            slon = [geoc.latitude]
        elif level in ['station']:
            for station in net:
                slat = station.latitude
                slon = station.longitude

    gcp2map(slat, slon, elat, elon, map=map, projection=proj)
    plt.show()
    return


def plot_inv_cat(invs=None, cats=None, axis=None, mark=None, seg=None,
                 pick=None, colorlist=None, catlabels=None, invlabels=None,
                 projection='moll', lon0=175.,
                 resolution='c'):

    if not invs and not cats:
        msg = 'Give an inventory or a catalog to plot'
        raise IOError(msg)

    if mark and type(mark) is not list:
        mark = [mark]
    else:
        mark = []

    m = Basemap(projection=projection, lon_0=lon0, resolution=resolution,
                ax=axis)
    m.drawmapboundary(fill_color='white')
    m.fillcontinents(color='#A9A9A9', lake_color='white')
    if resolution != 'c':
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(np.arange(-90., 120., 30.), labels=[0, 1, 0, 0],
                        linewidth=0.8)
    m.drawparallels(np.arange(-90., 120., 30.), linewidth=0.8)
    m.drawmeridians(np.arange(0., 420., 30.), linewidth=0.8)

    if invs:
        if type(invs) is not list:
            invs = [invs]
        for inv in invs:
            for net in inv:
                for s in net:
                    lon = s.longitude
                    lat = s.latitude
                    if s.code in mark:
                        m.scatter(lon, lat, 200, latlon=True,
                                  marker='v', color='yellow',
                                  zorder=11)

                    elif seg:
                        for p in seg:
                            if s.code in p.stats.station:
                                m.scatter(lon, lat, 200, latlon=True,
                                          marker='v', color='#e57248',
                                          zorder=11)
                        else:
                            m.scatter(lon, lat, 200, latlon=True,
                                      marker='v', color='#bc8f8f',
                                      zorder=10)
                    else:
                        m.scatter(lon, lat, 200, latlon=True,
                                  marker='v', color='#e57248',
                                  zorder=11)
                    # plt.annotate(s.code, xy=m(lon-2, lat),  xycoords='data',
                    #              zorder=11)
    if cats:
        if type(cats) is not list:
            cats = [cats]
        if colorlist is None:
            colormap = iter(cm.jet(np.linspace(0, 1, len(cats))))
        else:
            colormap = iter(colorlist)

        points = []
        i_old = -1
        for i, cat in enumerate(cats):
            c = next(colormap)
            for event in cat:
                lon = event.origins[0].longitude
                lat = event.origins[0].latitude
                _m = m.scatter(lon, lat, 80, latlon=True,
                               marker='o', color=c,
                               zorder=11, edgecolors='k')
                if i != i_old:
                    points.append(_m)
                i_old = i

        if catlabels is not None:
            plt.legend(points, catlabels, ncol=len(catlabels),
                       loc='upper center', bbox_to_anchor=(0.5, -0.05))
    return


def plot_point_on_map(lat, lon, projection='moll', background=None,
                      basemap=True, lat0=None, lon0=None, **plotargs):
    """
    Tool to plot a given latitude-longitude point on a map

    projection:
    https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html
    """
    if basemap is True:
        m = Basemap(projection='nsper', resolution=None,
                    lat_0=lat0, lon_0=lon0)
        m.shadedrelief()
        m.scatter(lon, lat, 800, latlon=True, zorder=500, **plotargs)
        return

    if 'ax' not in plotargs or 'fig' not in plotargs:
        fig, ax = plt.subplots()

    if projection in ['nsper', 'NearsidePerspective']:
        ax = plt.axes(projection=ccrs.NearsidePerspective())
    elif projection == 'ESRI':
        ax = plt.axes(projection=ShadedReliefESRI().crs)
        # ax.set_extent([-22, -15, 63, 65])
        ax.add_image(ShadedReliefESRI(), 8)
    else:
        ax = plt.axes(projection=ccrs.Mollweide())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()
    # ax.stock_img()
    if background is None:
        ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
        ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
        ax.add_feature(cartopy.feature.COASTLINE)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
        ax.add_feature(cartopy.feature.RIVERS)
        ax.coastlines()
        ax.gridlines()

    elif background is not None and projection != 'ESRI':
        ax.background_img(name=background, **plotargs)

    ax.plot(lon, lat, transform=ccrs.PlateCarree(),
            zorder=12, **plotargs)
    return fig, ax


def plot_inv_cat_c2py(invs=None, cats=None, axis=None, mark=None, seg=None,
                      colorlist=None):

    if not invs and not cats:
        msg = 'Give an inventory or a catalog to plot'
        raise IOError(msg)

    if mark and type(mark) is not list:
        mark = [mark]
    elif not mark:
        mark = []

    if axis:
        ax = axis
    else:
        ax = plt.axes(projection=ccrs.Mollweide())
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

    if invs:
        lat_used = []
        lon_used = []
        lat_unused = []
        lon_unused = []
        lon_all = []
        lat_all = []
        if type(invs) is not list:
            invs = [invs]
        for inv in invs:
            for net in inv:
                for s in net:
                    lon = s.longitude
                    lat = s.latitude
                    if s.code in mark:
                        ax.plot(lon, lat, transform=ccrs.PlateCarree(),
                                marker='v', color=marked_station_c,
                                zorder=12, mew=marker_size)

                    elif seg:
                        s_plotted = []
                        for p in seg:
                            if s.code in p.stats.station:
                                if s.code in s_plotted:
                                    continue
                                lat_used = np.append(lat_used, [lat, np.nan])
                                lon_used = np.append(lon_used, [lon, np.nan])
                                s_plotted.append(s.code)
                            else:
                                lat_unused = np.append(lat_unused,
                                                       [lat, np.nan])
                                lon_unused = np.append(lon_unused,
                                                       [lon, np.nan])
                    else:
                        lat_all = np.append(lat_all, [lat, np.nan])
                        lon_all = np.append(lon_all, [lon, np.nan])

        if len(lon_used) != 0:
            ax.plot(lon_used, lat_used, transform=ccrs.PlateCarree(),
                    marker='v', color=used_station_c,
                    zorder=11, mew=marker_size)
        if len(lon_unused) != 0:
            ax.plot(lon_unused, lat_unused, transform=ccrs.PlateCarree(),
                    marker='v', color=unused_station_c,
                    zorder=10, mew=marker_size)
        if len(lon_all) != 0:
            ax.plot(lon_all, lat_all, transform=ccrs.PlateCarree(),
                    marker='v', color=all_station_c,
                    zorder=11, mew=marker_size)

    if cats:
        if type(cats) is not list:
            cats = [cats]
        if colorlist is None:
            colormap = iter(cm.jet(np.linspace(0, 1, len(cats))))
        else:
            colormap = iter(colorlist)

        lon_all = []
        lat_all = []
        for cat in cats:
            c = next(colormap)
            for event in cat:
                lon = event.origins[0].longitude
                lat = event.origins[0].latitude
                lat_all = np.append(lat_all, [lat, np.nan])
                lon_all = np.append(lon_all, [lon, np.nan])

        ax.plot(lon, lat, transform=ccrs.PlateCarree(),
                marker='o', color=c,
                zorder=11, mew=marker_size)

    return ax


def gcp2map(slat, slon, elat, elon, axis=None, cmt=None,
            projection='moll', map=None, proximity=None):

    """
    Description follows
    """

    # Setting up the map projection
    if projection not in ['kav7', 'nsper', 'moll']:
        projection = 'moll'

    difflon = abs(slon - elon)
    if difflon > 180:
        lon0 = int((slon+elon)/2. + 180)
    else:
        lon0 = int((slon+elon)/2.)

    if projection == 'kav7' or projection == 'moll':
        if axis:
            m = Basemap(projection=projection, lon_0=lon0, resolution='c',
                        ax=axis)
        else:
            m = Basemap(projection=projection, lon_0=lon0, resolution='c')

    elif projection == 'nsper':
        lat0 = (slat+elat)/2.
        if axis:
            m = Basemap(projection='nsper', resolution=None,
                        lat_0=lat0, lon_0=lon0, ax=axis)
        else:
            m = Basemap(projection='nsper', resolution=None,
                        lat_0=lat0, lon_0=lon0)

    # Background of the map
    if map == 'bluemarble':
        m.bluemarble()
        m.drawgreatcircle(elon, elat, slon, slat, linewidth=3,
                          color='white', zorder=1)
    elif map == 'shadedrelief':
        m.shadedrelief()
        m.drawgreatcircle(elon, elat, slon, slat, linewidth=3,
                          color='black', zorder=1)
    elif map == 'etopo':
        m.etopo()
        m.drawgreatcircle(elon, elat, slon, slat, linewidth=3,
                          color='black', zorder=1)
    else:
        m.drawmapboundary(fill_color='lightgray')
        m.fillcontinents(color='white', lake_color='lightgray')
        m.drawcoastlines()
        #m.drawrivers()
        m.drawcountries(linestyle=':')

        m.drawparallels(np.arange(-90., 120., 30.), color='gray')
        if projection == 'moll':
            m.drawmeridians(np.arange(0., 420., 60.), color='gray')
        else:
            m.drawmeridians(np.arange(0., 420., 60.), color='gray')
        m.drawgreatcircle(elon, elat, slon, slat, linewidth=3,
                          color='red', zorder=1)

    # Plotting station

    if proximity:
        d_km = gps2dist_azimuth(slat, slon, elat, elon)[0]/1000.
        d_deg = kilometer2degrees(d_km)
        if d_deg <= proximity:
            scolor = 'red'
        else:
            scolor = 'yellow'
    else:
        scolor = 'yellow'

    m.scatter(slon, slat, 500, latlon=True, marker='v', color=scolor,
              edgecolors= "black", zorder=10)

    # Plotting event
    if cmt:
        axins = inset_axes(axis, .8, .8)
        xmax = axis.get_xlim()[1]
        ymax = axis.get_ylim()[1]
        xpt, ypt = m(elon, elat)
        h = .17
        w = .075
        xpt = xpt/xmax - w/2.
        ypt = ypt/ymax - h/2.
        ip = InsetPosition(axis, [xpt, ypt, w, h])
        axins.set_axes_locator(ip)
        beachball(cmt, linewidth=0.5, axis=axins)

    else:
        m.scatter(elon, elat, 300, latlon=True, marker='*', color='#CC0000',
                  zorder=15)
    return axis


def gcpmap_c2py(slat, slon, elat, elon, axis=None, cmt=None,
                projection='moll', proximity=None):

    """
    Description follows
    """
    if cmt is not None:
        cmt = cmt[0:6]
    if axis:
        ax = axis
    else:
        difflon = abs(slon - elon)
        if difflon > 180:
            lon0 = int((slon+elon)/2. + 180)
        else:
            lon0 = int((slon+elon)/2.)
        ax = plt.axes(projection=ccrs.Mollweide(central_longitude=lon0))
        # Plotting station
        ax.set_global()
        ax.add_feature(cartopy.feature.LAND, facecolor=land_color)
        ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
        ax.add_feature(cartopy.feature.COASTLINE)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
        ax.add_feature(cartopy.feature.RIVERS)
        ax.coastlines()
        ax.gridlines()
    if proximity:
        d_km = gps2dist_azimuth(slat, slon, elat, elon)[0]/1000.
        d_deg = kilometer2degrees(d_km)
        if d_deg <= proximity:
            scolor = 'red'
        else:
            scolor = 'yellow'
    else:
        scolor = 'yellow'

    ax.plot(slon, slat, transform=ccrs.PlateCarree(), marker='v', color=scolor,
            zorder=10, mew=marker_size)

    # Currently plotted as part of the main plot, [0.88, 0.05, 0.1, 0.1] refers
    # to the the coordinate system of the main figure

    # The Beachball is slightly too big, needs fixing
    # sub_ax = plt.axes([0.55, 0.23, 0.1, 0.1], projection=ccrs.PlateCarree())
    # sub_ax.outline_patch.set_alpha(0.0)
    # sub_ax.background_patch.set_fill(False)
    # if cmt:
    #     beachball(cmt, axis=sub_ax)

    # else:
    ax.plot(elon, elat, transform=ccrs.PlateCarree(), marker='*',
            color='#0000FF', zorder=15, mew=marker_size)
    ax.plot([slon, elon], [slat, elat], transform=ccrs.Geodetic(),
            color='#CC0000')
    return ax


def plot_cst_map(infile, outfile=None):

    with open(infile, 'r') as fh:
        c_lat_lon = np.genfromtxt(fh)

    projection = 'moll'
    lon0 = 0.
    m = Basemap(projection=projection, lon_0=lon0, resolution='c')

    c = c_lat_lon.transpose()
    lats = c[0]
    lons = c[1]
    data = c[2]
    x, y = m(lons, lats)
    m.scatter(x, y, data, zorder=15)
    m.drawmapboundary(fill_color='white')
    m.fillcontinents(color='#dfdfdf', lake_color='white')
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(-90., 120., 30.), labels=[0, 1, 0, 0])
    m.drawmeridians(np.arange(0., 420., 60.))

    return
