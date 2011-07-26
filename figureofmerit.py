# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import polarizer2 as pol


mirrordef = {
        'names': ['B4C', 'Mo'],
        'thickness': [3, 100],
        'energy': 60,
        }
mirror = pol.mirror(mirrordef)

fig = plt.figure()
# transmission, figure of merit
ax1 = fig.add_subplot(111)
# degree of circ. pol
ax2 = ax1.twinx()
angles = np.arange(np.pi / 4, np.pi / 2, .005)
R_s, R_p, delta = np.array([mirror.get_Rd(a) for a in angles]).transpose()

beta = np.arctan(R_p**4 / R_s**4)
R = np.sqrt((np.cos(beta) * R_p ** 4)**2 + (np.sin(beta) * R_s ** 4)**2)
DoP = np.sin(4 * delta)
FoM = R * DoP ** 2
# print np.sin(beta) * R_s ** 4 - np.cos(beta) * R_p ** 4
# print np.rad2deg(beta)

# angle-direction is reversed (grazing incidence vs. normal incidence)
pangle = np.rad2deg(np.pi / 2 - angles)
pl1 = ax2.plot(pangle, DoP, lw=2, color='red')

pl2 = ax1.plot(pangle, FoM, lw=2, color='indigo')
Rplot = ax1.plot(pangle, R, lw=2, color='darkgreen')
# Rsplot = ax1.plot(pangle, R_s**4, lw=2, color='blue')
# Rpplot = ax1.plot(pangle, R_p**4, lw=2, color='green')

ax1.set_xlabel('incidence angle (deg)')
ax2.set_ylabel('degree of circular polarization', color='darkred')
ax1.set_ylabel('total transmission', color='darkgreen')

# annotations
arrowprops = {'arrowstyle': 'simple', 'color': 'black'}
target = np.where(FoM == FoM.max())
m = ax2.axvline(pangle[target], color='black', ls='-.')
DoPtext = '$P_c=%.2f$' % DoP[target]
DoPxy = (pangle[target], DoP[target])
DoPtxy = (pangle[target] + 2.5, DoP[target] - .1)
ax2.annotate(DoPtext, DoPxy, DoPtxy, arrowprops=arrowprops)

Rtext = '$T=%.2f$' % R[target]
Rxy = (pangle[target], R[target])
Rtxy = (pangle[target] + 2.5, R[target] + .02)
ax1.annotate(Rtext, Rxy, Rtxy, arrowprops=arrowprops)

# DoPmax = np.where(DoP == DoP.max())
# target90 = np.where(DoP[:DoPmax] == min(DoP[:DoPmax], key=lambda x: abs(x-.89)))
# ax2.axvline(pangle[target90], color='black', ls='-.')

leg = plt.legend((pl1, Rplot, pl2), ('$P_c$', '$T$', u'$T·P_c^2$'), 
        loc='upper right', fancybox=True,)
leg.legendPatch.set_alpha(0.7)

# fig.savefig('figureofmerit.pdf')
fig.show()
# vim: foldmethod=marker
