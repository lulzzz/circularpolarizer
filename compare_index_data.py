# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
import polarizer2 as pol


def plotmirror(ind):
    # substrate thickness is not needed, but set to -1 to avoid IndexError
    mirrordef = {
            'names': ['B4C', ind[0], 'SiO2'],
            'thickness': [3, 50, -1],
            'energy': 60,
            }
    mirror= pol.mirror(mirrordef)
    Rd = np.array([mirror.get_Rd(a) for a in angles])
    beta = np.arctan(Rd[:,0] / Rd[:,1])
    R = (np.sin(beta) * Rd[:,0]) ** 4
    DoP = np.sin(4 * Rd[:,2])
    FoM = R * DoP ** 2

    ax.plot(pangle, DoP, lw=2, color='red', ls=ind[1])
    ax2.plot(pangle, FoM, lw=2, color='indigo', ls=ind[1])
    ax2.plot(pangle, R, lw=2, color='darkgreen', ls=ind[1])


fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.set_xlabel('incidence angle (deg)')
ax.set_ylabel('degree of circular polarization', color='darkred')
ax2.set_ylabel('total transmission', color='darkgreen')

angles = np.arange(np.pi / 4, np.pi / 2, .005)
# angle-direction is reversed (grazing incidence vs. normal incidence)
pangle = np.rad2deg(np.pi / 2 - angles)

indexfiles = (('Mo', '-'), ('Mo_palik', '--'), ('Mo_llnl_cxro', '-.'))
[plotmirror(ind) for ind in indexfiles]

ax.set_xlim(0, 45)
ax.set_ylim(0, 1)
ax2.set_ylim(0, .25)

# fake lines for the legend:
red = plt.plot(2, 2, color='red', lw=2)
darkgreen = plt.plot(2, 2, color='darkgreen', lw=2)
indigo = plt.plot(2, 2, color='indigo', lw=2)
cxro = plt.plot(2, 2, color='black', lw=2, ls='-')
llnl = plt.plot(2, 2, color='black', lw=2, ls='-.')
palik = plt.plot(2, 2, color='black', lw=2, ls='--')
leg = ax.legend((red, darkgreen, indigo), ('$P_c$', '$T$', u'$T·P_c^2$'), 
        loc='upper right', fancybox=True,)
leg2 = ax2.legend((cxro, llnl, palik), ('$CXRO$', '$LLNL$', '$Palik$'),
        loc='lower right', fancybox=True,)

leg.legendPatch.set_alpha(0.7)
leg2.legendPatch.set_alpha(0.7)

fig.savefig('compare_index_data.pdf')
# fig.show()
# vim: foldmethod=marker
