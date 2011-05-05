# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, hbar
from os import getcwd

# {{{ definitions

# Müller-Matrix
M = lambda phi,delta,rp,rs: (rp**2+rs**2)/2*np.matrix([
    [1,-np.cos(2*phi),0,0], \
    [-np.cos(2*phi),1,0,0], \
    [0,0,np.sin(2*phi)*np.cos(delta),np.sin(2*phi)*np.sin(delta)], \
    [0,0,-np.sin(2*phi)*np.sin(delta),np.sin(2*phi)*np.cos(delta)]])

# Stokes-Vektor
S = lambda alpha: np.array([1, np.cos(2*alpha), np.sin(2*alpha), 0])

# Kramers-Kronig-Relation: 
# calculate the imaginary part of reflectivity (=phaseshift) at given circ.frequency
# omega: given frequency
# _omega: list of all frequencies to be 'integrated' (=sum) over
# _refl: corresponding list of reflectivities
# kk_phase = lambda omega, _omega, _refl: \
        # omega/np.pi * sum([ np.log(_refl[i]/_refl[np.where(_omega == omega)])/(om**2 - omega**2) \
        # for i, om in enumerate(_omega) if om!=omega ])

# convert Energy(eV) -> ω
eV2omega = lambda eV: eV*e/hbar

# import dat-files from CXRO-bilayer-refl-calculator, returns angles, reflectivities
def importindex(filename):        
    return np.genfromtxt(filename, skip_header=2, unpack=True)

# put it all together...
# def phaseshift(filename, energy):
    # _energy, reflectivity = importdat(filename)
    # phaseshift = kk_phase(eV2omega(energy), eV2omega(_energy), reflectivity)
    # return phaseshift

# refl_i = complex reflection coeff. ~ r*exp(iφ)
refl_s = lambda eps, th: eps*np.sin(th) - np.sqrt((eps-np.cos(th)**2))\
        /(eps*np.sin(th) + np.sqrt(eps-np.cos(th)**2))

refl_p = lambda eps, th: np.sin(th) - np.sqrt((eps-np.cos(th)**2))\
        /(np.sin(th) + np.sqrt(eps-np.cos(th)**2))

# get phase φ from complex number C with abs(C)=r: C=r*exp(іφ)
cplx_phase = lambda cplx: np.arctan(cplx.real/cplx.imag)

# convert complex refractive index N=n+ik to complex dielectric constant
n2eps = lambda n, k: np.complex(n**2 - k**2, 2*n*k)

# closest match-finder
closest = lambda liste, match: min(list(liste), key = lambda i: abs(i-match))
# }}}

startdir = '/home/dscran/Documents/promotion/circularpolarizer'
# phasediff = [ phaseshift('%s/data/deg%02d_p.dat' % (startdir, angle), energy) \
        # - phaseshift('%s/data/deg%02d_s.dat' % (startdir, angle), energy) for angle in angles ]

target_energy = 60              #eV

# the refractive index data files contain [E(eV), δ, β] with n=1-δ-iβ
energy, delta, beta = importindex(startdir + '/data/index_B4C.dat')

delta = 1 - delta[ energy == closest(energy, target_energy) ]
beta = -beta[ energy == closest(energy, target_energy) ]

eps = n2eps( delta, beta )

theta = np.arange(0, np.pi/2, .01)

phi_s = np.array( [ cplx_phase( refl_s(eps, th) ) for th in theta ] )
phi_p = np.array( [ cplx_phase( refl_p(eps, th) ) for th in theta ] )

D_phi = phi_s - phi_p

# {{{ Plotten

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(np.rad2deg(theta), np.rad2deg(D_phi), 'g-')
# plt.plot(angles, reflectivities_s, 'b-')
fig.savefig(startdir+'/phase.pdf')

# }}}

# vim: foldmethod=marker
