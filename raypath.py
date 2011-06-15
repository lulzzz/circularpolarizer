# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


pi = np.pi
# capital letters denote points, small ones stand for direction of straights.
class straightline(object):
    def __init__(self, basepoint, angle):
        self.bp = basepoint
        self.angle = angle
        self.d = np.array([np.cos(angle), np.sin(angle)])

    def goto(self, l):
        return self.bp + l * self.d

    # intersection with other straightline
    def isct(self, line1):
        E0x, E0y = self.bp
        e0x, e0y = self.d
        E1x, E1y = line1.bp
        e1x, e1y = line1.d
        # solution of l for straight==line1
        l = -(e1x * E1y - e1y * E1x + E0x * e1y - E0y * e1x) /\
                (e0x * e1y - e0y * e1x)
        return l

    # reflection on a mirror straight, returns basepoint and angle of
    # reflected straightline
    def refl(self, mirror):
        return self.goto(self.isct(mirror)), 2 * mirror.angle - self.angle

    def mirrorpath(self, mirrorlist):
        ray = self
        raypath = [self.bp]
        for m in mirrorlist:
            bp, angle = ray.refl(m)
            ray = straightline(bp, angle)
            raypath.append(bp)
        raypath.append(ray.goto(mirrorlist[-1].bp[0]+150-bp[0]))
        return np.array(raypath)

# close existing figure, if any
try:
    plt.close(plt.gcf())
except:
    pass

# angle of incidence (mirrors) and divergence angle of ray
a = np.deg2rad(22.5)
omega = np.deg2rad(3)
# surface-normal distance of mirror pairs (mm)
d = 50
# distance between the two mirror-pairs (mm)
l = 300
# mirrorsize
ms = 120
# focus point = origin of light-rays
f = [-100, 0]

def genmirrors(a):
    return [
            straightline([0, 0], a),
            straightline([d * np.cos(pi / 2 + a), d * np.sin(pi / 2 + a)], a),
            straightline([l + d * np.sin(a), d * np.cos(a)], -a),
            straightline([l, 0], -a)
            ]

fig = plt.figure()
pic = fig.add_subplot(111)
fig.subplots_adjust(bottom=.3)
pic.axis('equal')
mirrors = genmirrors(a)
pathplots, mirrorplots = [], []
for o in [-omega, 0, omega]:
    path = straightline(f, o).mirrorpath(mirrors).T
    c = 'red' if o else 'darkgreen'
    pathplots.append(pic.plot(path[0], path[1], '-', color=c, lw=2)[0])
for m in mirrors:
    M = np.array([m.goto(-ms), m.goto(ms)])
    mirrorplots.append(pic.plot(M[:,0], M[:,1], '-', color='black', lw=2)[0])

axcolor = 'lightgoldenrodyellow'
axangle = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axomega  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
sangle = Slider(axangle, 'incidence', 0, 45.0, valinit=22.5)
somega = Slider(axomega, 'divergence', 1, 10.0, valinit=1)

def update(val):
    omega = np.deg2rad(somega.val)
    a = np.deg2rad(sangle.val)
    mirrors = genmirrors(a)
    for n,o in enumerate([-omega, 0, omega]):
        path = straightline(f, o).mirrorpath(mirrors).T
        pathplots[n].set_data(path[0], path[1])
    for n,m in enumerate(mirrors):
        M = np.array([m.goto(-ms), m.goto(ms)])
        mirrorplots[n].set_data(M[:,0], M[:,1])
        fig.show()
sangle.on_changed(update)
somega.on_changed(update)
fig.show()
