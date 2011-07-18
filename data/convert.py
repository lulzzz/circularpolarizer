# -*- coding: utf-8 -*-

# simply convert 'foreign' refractive-indices-data to the format expected by
# 'polarizer2.py'

# Mo_llnl_cxro MO_windt88 Mo_palik
# input format is:     λ(Å)    n    k
# output:              E(eV)   δ    β       with n=1-δ; k=β

import numpy as np
from scipy.constants import h, c, e


def convert(i):
    A, n, k = np.genfromtxt(i + '.nk', comments=';', unpack=True)
    E = angstrom2eV(A)
    o = np.array((E, 1 - n, k)).transpose()
    np.savetxt('n_%s.dat' % i, o)


angstrom2eV = lambda A: h * c / (e * 1e-10 * A)


infiles = ['Mo_llnl_cxro', 'Mo_windt88', 'Mo_palik']
[convert(i) for i in infiles]
