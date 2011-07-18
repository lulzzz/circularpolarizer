# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import polarizer2 as pol


# substrate thickness is not needed, but set to -1 to avoid IndexError
mirrordef = {
        'names': ['B4C', 'Mo', 'SiO2'],
        'thickness': [3, 50, -1],
        'energy': 60,
        }
mirror = pol.mirror(mirrordef)

fig = plt.figure()
ax = fig.add_subplot(111)
# angles = np.arange(np.pi / 4, np.pi / 2, .005)
angles = np.array([np.deg2rad(65)])
Rd = np.array([mirror.get_Rd(a) for a in angles])

beta = np.arctan(Rd[:,0] / Rd[:,1])
R = (np.sin(beta) * Rd[:,0]) ** 4
DoP = np.sin(4 * Rd[:,2])
FoM = R * DoP ** 2

# angle-direction is reversed (grazing incidence vs. normal incidence)
pangle = np.rad2deg(np.pi / 2 - angles)
pl1 = ax.plot(pangle, DoP, lw=2, color='red')
m = ax.axvline(pangle[FoM == FoM.max()], color='black', ls='-.')

ax2 = ax.twinx()
pl2 = ax2.plot(pangle, FoM, lw=2, color='indigo')
Rplot = ax2.plot(pangle, R, lw=2, color='darkgreen')

# ax.set_ylim(ymax=190)
ax.set_xlabel('incidence angle (deg)')
ax.set_ylabel('degree of circular polarization', color='darkred')
ax2.set_ylabel('total transmission', color='darkgreen')

# annotations
arrowprops = {'arrowstyle': 'simple', 'color': 'black'}
target = np.where(FoM == FoM.max())
DoPtext = '$P_c=%.2f$' % DoP[target]
DoPxy = (pangle[target], DoP[target])
DoPtxy = (pangle[target] + 4, DoP[target] - .1)
ax.annotate(DoPtext, DoPxy, DoPtxy, arrowprops=arrowprops)

Rtext = '$T=%.2f$' % R[target]
Rxy = (pangle[target], R[target])
Rtxy = (pangle[target] + 4, R[target] + .02)
ax2.annotate(Rtext, Rxy, Rtxy, arrowprops=arrowprops)

leg = plt.legend((pl1, Rplot, pl2), ('$P_c$', '$T$', u'$T·P_c^2$'), 
        loc='upper right', fancybox=True,)
leg.legendPatch.set_alpha(0.7)

# fig.savefig('figureofmerit.pdf')
# fig.show()
# vim: foldmethod=marker
