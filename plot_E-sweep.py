# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np


def makename(names):
    fname = ''
    for s in names:
        fname += s + '-'
    return fname.rstrip('-')


# change the fonts
f = 'STIXGeneral'
font = {'family': f}
mfont = {'default': 'regular',
'fontset': 'stixsans'}
plt.rc('font', **font)
plt.rc('mathtext', **mfont)

mnames = [['B4C', 'Mo'], ['Rh'], ['Ru'], ['TiO', 'Ti'], ['W'], ['Pd']]
mnames = [['B4C', 'Mo'], ['Ru'], ['TiO', 'Ti']]

# close existing figure, if any
try:
    plt.close(plt.gcf())
except:
    pass

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for i, name in enumerate(mnames):
    fname = makename(name)
    x, R_s, R_p, angle = np.genfromtxt('data/interpolated_E-sweep_' + fname + '.txt', unpack=True)
    R_s_ges = R_s ** 4
    R_p_ges = R_p ** 4
    beta = np.arctan(R_p_ges / R_s_ges)
    ax1.plot(x, R_s_ges * np.sin(beta), 'o', label=u'%s' % fname)
    # ax1.plot(x, R_p_ges * np.cos(beta), 'o', ms=3)
    ax2.plot(x, np.rad2deg(beta), 'o', label=u'%s' % fname)
    # ax2.plot(x, 90 - np.rad2deg(angle), 'o', label=u'%s' % fname)

ax1.legend(loc='upper left', numpoints=1)
# ax2.legend(loc='upper left', numpoints=1)
ax2.set_ylim(0)
ax1.set_xlabel('energy (eV)')
ax2.set_xlabel('energy (eV)')
ax1.set_ylabel(u'reflectivity')
ax2.set_ylabel(u'beta (deg)')
# ax2.set_ylabel(u'required grazing angle of incidence (deg)')
fig.suptitle(u'reflectivity after 4 mirrors for total phaseshift of 90°')

# fig.savefig('E-sweep.pdf')
fig.show()
