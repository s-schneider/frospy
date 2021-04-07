import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import math

pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt


# Displacement functions
def horizontal_displ(npts_theta=36, npts_phi=9, lat='N', displ=0):
    """
    Calculates horizontal displacement by rotating all points on
    latitude lat by the value displ (adding it to the range of 0 to 2pi).

    npts_theta : Number of longitude grid points
    npts_phi   : Number of Latitude grid points
    lat        : Latitude that will be rotated, 'N', 'S', 'both'
                 or value in range of 0 to 2pi
    displ      : amount of displacement

    """
    r = 5
    _theta = np.linspace(0 + displ, 2 * pi + displ, npts_theta)

    if lat == 'N':
        _phi = np.linspace(0, pi / 2, npts_phi, endpoint=False)
    elif lat == 'S':
        _phi = np.linspace(pi / 2, pi, npts_phi)
    elif lat == 'both':
        _phi = np.linspace(0, pi, npts_phi)
    else:
        _phi = lat

    phi, theta = np.meshgrid(_phi, _theta)

    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    return x, y, z


def vertical_displ(npts_theta=36, npts_phi=9, lat=0, displ=0):
    """
    Calculates vertical displacement by changing the radius of all points on
    latitude lat by the value displ.

    npts_theta : Number of longitude grid points
    npts_phi   : Number of Latitude grid points
    lat        : value in range of 0 to 2pi
    displ      : amount of displacement

    """
    r = 5
    _theta = np.linspace(0, 2 * pi, npts_theta)
    _phi = lat
    phi, theta = np.meshgrid(_phi, _theta)

    x = (r + displ) * sin(phi) * cos(theta)
    y = (r + displ) * sin(phi) * sin(theta)
    z = (r + displ) * cos(phi)

    return x, y, z


def generate_grid_lines(lon, npts_theta=36, npts_phi=9):
    """
    Generates gridlines for a given longitude
    """
    r = 5
    _theta = lon
    _phi = np.linspace(0, 2 * pi, npts_phi)
    phi, theta = np.meshgrid(_phi, _theta)

    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    return x, y, z


# Distortion
def Tmovement(t, pattern, scatters, npts_theta, npts_phi):
    """
    Movement funciton for animation for toroidal mode pattern
    """
    # Distort the sphere here with pattern
    for i, (lat, scatter) in enumerate(scatters[0].items()):
        _t = t * pattern[i]

        XYZ = horizontal_displ(npts_theta=npts_theta, npts_phi=npts_phi,
                               lat=lat, t=_t)
        xyz = np.vstack([XYZ[0].ravel(), XYZ[1].ravel(), XYZ[2].ravel()])
        # Update points on the sphere
        for i, axis in enumerate(scatter._offsets3d):
            offset = axis
            for j, point in enumerate(xyz[i]):
                offset.data[j] = point
            scatter._offsets3d[i]
    plt.draw()
    return scatters


def Smovement(t, pattern, scatters, npts_theta, npts_phi):
    """
    Movement funciton for animation for spheroidal mode pattern
    """
    # Distort the sphere here with pattern
    for i, (lat, scatter) in enumerate(scatters[0].items()):
        _t = t * pattern[i]

        XYZ = vertical_displ(npts_theta=npts_theta,
                             npts_phi=npts_phi,
                             lat=lat, displ=_t)
        xyz = np.vstack([XYZ[0].ravel(), XYZ[1].ravel(), XYZ[2].ravel()])
        # Update points on the sphere
        for i, axis in enumerate(scatter._offsets3d):
            offset = axis
            for j, point in enumerate(xyz[i]):
                offset.data[j] = point
            scatter._offsets3d[i]
    plt.draw()
    return scatters


def STmovement(t, pattern, scatters, npts_theta, npts_phi):
    """
    Movement funciton for animation for superposition of
    toroidal and spheroidal mode pattern
    """
    # Distort the sphere here with pattern
    for i, (lat, scatter) in enumerate(scatters[0].items()):
        _T = t * pattern[0][i]
        _S = t * pattern[1][i]

        TXYZ = horizontal_displ(npts_theta=npts_theta, npts_phi=npts_phi,
                                lat=lat, displ=_T)
        SXYZ = vertical_displ(npts_theta=npts_theta,
                              npts_phi=npts_phi,
                              lat=lat, displ=_S)
        xyz = np.vstack([(SXYZ[0].ravel() + TXYZ[0].ravel()) / 2.,
                         (SXYZ[1].ravel() + TXYZ[1].ravel()) / 2.,
                         (SXYZ[2].ravel() + TXYZ[2].ravel()) / 2.])
        # Update points on the sphere
        for i, axis in enumerate(scatter._offsets3d):
            offset = axis
            for j, point in enumerate(xyz[i]):
                offset.data[j] = point
            scatter._offsets3d[i]
    plt.draw()
    return scatters


def create_plot(frames, interval=100, npts_theta=27, npts_phi=13,
                surface=False):
    """
    Creates the figure with a 3D scatter plot, hemispheres are
    N: blue
    S: orange
    equator: grey

    returns: scatters, fig and ax
    """
    fig = plt.figure()
    fig.set_size_inches(9, 9)
    ax = fig.add_subplot(111, projection='3d', label='axes1')
    ax.view_init(18, 20)
    # Parametric eq for a distorted globe (for demo purposes)
    scatters = [{}]
    Nlats = np.linspace(0, pi / 2, math.floor(npts_phi / 2.), endpoint=False)
    Slats = np.linspace(pi / 2, pi, math.ceil(npts_phi / 2.))
    lats = np.hstack((Nlats, Slats))
    for lat in lats:
        xyz = horizontal_displ(npts_theta=npts_theta, npts_phi=npts_phi,
                               lat=lat, displ=0)
        if lat > pi / 2:
            scatters[0][lat] = ax.scatter(xyz[0], xyz[1], xyz[2], c='#ff7f0e')
        elif lat == pi / 2:
            scatters[0][lat] = ax.scatter(xyz[0], xyz[1], xyz[2], c='#808080')
        else:
            scatters[0][lat] = ax.scatter(xyz[0], xyz[1], xyz[2], c='#1f77b4')
        # plot meridians
        # ax.plot(xyz[0].T[0], xyz[1].T[0], zs=xyz[2].T[0], c='#808080')

    # if surface is True:
    #     xyz = generate_grid(npts_theta=npts_theta, npts_phi=npts_phi,
    #                         hemisphere='both')
    #     surf = ax.plot_surface(xyz[0], xyz[1], xyz[2])
    #     return scatters, [surf], fig

    # lons = np.linspace(0, 2 * pi, npts_theta)
    # for lon in lons:
    #     xyz = generate_grid_lines(lon=lon, npts_theta=npts_theta,
    #                               npts_phi=npts_phi*10)
    #     ax.plot(xyz[0][0], xyz[1][0], zs=xyz[2][0], c='#808080')
    return scatters, fig, ax


def init_plot(mode='T', grid='off', axis='off',
              npts_theta=27, npts_phi=13, interval=100):
    """
    Initialize the plot and movement functions for the animation.
    """
    _x = np.linspace(0, 2 * pi, 90)
    frames = np.sin(_x)
    scatters, fig, ax = create_plot(frames, interval, npts_theta, npts_phi)

    _x = np.linspace(0, 2 * pi, len(scatters[0]))
    if mode == 'T':
        pattern = np.sin(_x)
        movement = Tmovement
    elif mode == 'S':
        pattern = 0.8 * np.cos(_x)
        movement = Smovement
    else:
        pattern = [np.sin(_x), 0.8 * np.cos(_x)]
        movement = STmovement

    if grid == 'off':
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    if axis == 'off':
        ax.set_axis_off()
    return scatters, pattern, movement, frames, fig, ax


def main(save=True, interval=100, npts_theta=27, npts_phi=13,
         mode='T', grid='off', axis='off'):
    """
    Animation funtion
    mode       : 'T', 'S' or 'ST'
    npts_theta : Number of longitude grid points
    npts_phi   : Number of latitude grid points
    interval   : inverval in ms between each frame
    """
    out = init_plot(mode=mode, grid=grid,
                    axis=axis,
                    npts_theta=npts_theta,
                    npts_phi=npts_phi,
                    interval=interval
                    )
    scatters, pattern, movement, frames, fig, ax = out[:]
    ani = animation.FuncAnimation(fig, movement, frames=frames,
                                  fargs=(pattern, scatters,
                                         npts_theta, npts_phi),
                                  interval=interval, blit=False)

    if save:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800,
                        extra_args=['-vcodec', 'libx264'])
        ani.save('0{}2.mp4'.format(mode), writer=writer)


plt.ion()
main(mode='ST', save=True)
plt.show()
