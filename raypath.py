# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy import interpolate
import polarizer2 as pol

# {{{ functions...
pi = np.pi


def genmirrors(geometry):
    a = geometry['angle']
    d = geometry['normal']
    l = geometry['length']
    return [
            pol.straightline([0, 0], a),
            # TODO: french model: mirror shift fixed @ 37.5mm
            # mirrors at central-beam position
            pol.straightline([d * np.cos(pi / 2 + a),
                d * np.sin(pi / 2 + a)], a),
            pol.straightline([l + d * np.sin(a),
                d * np.cos(a)], -a),
            pol.straightline([l, 0], -a)
            ]


# align mirrors along a given raypath
def genmirrors2(ray, a):
    return [
            pol.straightline(ray[0], a),
            pol.straightline(ray[1], a),
            pol.straightline(ray[2], -a),
            pol.straightline(ray[3], -a)
            ]


# update plot when sliders are changed
def update(val):
    omega = np.deg2rad(somega.val)
    geometry['angle'] = np.deg2rad(sangle.val)
    mirrors = genmirrors(geometry)
    for n, o in enumerate([-omega, omega, 0]):
        path = pol.straightline(geometry['focus'], o).mirrorpath(mirrors).T
        pathplots[n].set_data(path[0], path[1])
        if not o:
            mirrors2 = genmirrors2(path.T[1:5], geometry['angle'])
    for n, m in enumerate(mirrors2):
        M = np.array([m.goto(-geometry['msize'][n]),
            m.goto(geometry['msize'][n])])
        mirrorplots[n].set_data(M[:, 0], M[:, 1])
        fig.show()
    m1.set_xdata(np.rad2deg(geometry['angle']))
    m2.set_xdata(np.rad2deg(geometry['angle'] - omega))
    m3.set_xdata(np.rad2deg(geometry['angle'] + omega))
    geotext.set_text('R_s = %.2f\nphase = %.2f'\
            % (R_s_fit(np.pi / 2 - geometry['angle']) ** 4,
                4 * np.rad2deg(phase_fit(np.pi / 2 - geometry['angle']))))

# }}}

# close existing figure, if any
try:
    plt.close(plt.gcf())
except:
    pass

# divergence angle of ray
omega = np.deg2rad(2)
# mirror settings: normal surface to surface distance of mirror pairs,
# longitudinal spacing, mirrorsize, focusposition, angle of incidence
geometry = {
        'normal': 30,
        'length': 230,
        'msize': (50, 50, 50, 50),
        'focus': (-100, 0),
        'angle': np.deg2rad(21)
        }

# substrate thickness is not needed, but set to -1 to get right visualization
mirrordef = {
        'names': ('B4C', 'Mo', 'SiO2'),
        'thickness': (3, 50, -1),
        'energy': 60,
        }

# {{{ plot
# {{{ beam visualization
fig = plt.figure()
pic = fig.add_subplot(211)
fig.subplots_adjust(top=.95, left=.07, bottom=.15, right=.75)
pic.axis('equal')
pic.set_xlabel('x (mm)')
pic.set_ylabel('y (mm)')
pic.set_xlim(-110, 360)
mirrors = genmirrors(geometry)
pathplots, mirrorplots = [], []
for o in [-omega, omega, 0]:
    path = pol.straightline(geometry['focus'], o).mirrorpath(mirrors).T
    c = 'red' if o else 'darkgreen'
    pathplots.append(pic.plot(path[0], path[1], '-', color=c, lw=2)[0])
    if not o:
        mirrors2 = genmirrors2(path.T[1:5], geometry['angle'])
for i, m in enumerate(mirrors2):
    M = np.array([m.goto(-geometry['msize'][i]), m.goto(geometry['msize'][i])])
    mirrorplots.append(pic.plot(M[:, 0], M[:, 1], '-', color='black', lw=2)[0])
# }}}

# {{{ polarization and reflectivity
mirror = pol.mirror(mirrordef)
angles = np.arange(np.pi / 4, np.pi / 2, .01)
Rd = np.array([mirror.get_Rd(a) for a in angles])
# probably overkill to interpolate those...maybe replace with simple
# closest-match finder
R_s_fit = interpolate.UnivariateSpline(angles, Rd[:, 0], s=0)
phase_fit = interpolate.UnivariateSpline(angles, Rd[:, 2], s=0)
pic2 = fig.add_subplot(212)
# angle-direction is reversed (grazing incidence vs. normal incidence)
delta_p = pic2.plot(np.rad2deg(np.pi / 2 - angles), 4 * np.rad2deg(Rd[:, 2]),
        lw=2, color='darkviolet')
pic2.set_ylim(ymax=190)
pic2.axhline(90, ls='-.', color='black')
pic2.set_xlabel('incidence angle (deg)')
pic2.set_ylabel('total phasediff (deg)')
m1 = pic2.axvline(np.rad2deg(geometry['angle']), color='darkgreen', lw=2)
m2 = pic2.axvline(np.rad2deg(geometry['angle'] - omega), color='red', lw=2)
m3 = pic2.axvline(np.rad2deg(geometry['angle'] + omega), color='red', lw=2)
pic3 = pic2.twinx()
pic3.set_ylabel('total reflectivity')
R = pic3.plot(np.rad2deg(np.pi / 2 - angles), Rd[:, 0] ** 4,
        lw=2, color='blue')
leg = plt.legend((delta_p, R), ('phase', 'reflectivity'), loc='lower right')
# }}}

# {{{ info window
info = plt.axes([.8, .15, .15, .80], axisbg='#bbbbbb')
info.set_xticks([])
info.set_yticks([])
info.text(.05, .30, mirror.info())
geotext = info.text(.05, .15, 'R_s = %.2f\nphase = %.2f'\
        % (R_s_fit(np.pi / 2 - geometry['angle']) ** 4,
            4 * np.rad2deg(phase_fit(np.pi / 2 - geometry['angle']))))
axcolor = 'white'
# sliders are (interactive) axes...
axangle = plt.axes([0.10, 0.02, 0.30, 0.03], axisbg=axcolor)
axomega = plt.axes([0.60, 0.02, 0.30, 0.03], axisbg=axcolor)
sangle = Slider(axangle, 'incidence', 0, 45.0,
        valinit=np.rad2deg(geometry['angle']))
somega = Slider(axomega, 'divergence', 0, 10.0,
        valinit=np.rad2deg(omega))
sangle.on_changed(update)
somega.on_changed(update)

fig.show()
# vim: foldmethod=marker
