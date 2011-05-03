# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

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
# calculate the imaginary part of reflectivity (=phaseshift) at angle phi
kk_phase = lambda phi, angles, refl: \
        -2/np.pi*sum([phi*refl[i]/(angle**2-phi**2) for i,angle in enumerate(angles) ])

# import dat-files from CXRO-bilayer-refl-calculator, returns angles, reflectivities
def importdat(filename):        
    return np.genfromtxt(filename, skip_header=3, usecols=[0,1], unpack=True)

# }}}

angles, reflectivities_p = importdat('xray_p.dat')
phases_p = np.array([kk_phase(np.deg2rad(phi), angles, reflectivities_p) for phi in angles])

angles, reflectivities_s = importdat('xray_s.dat')
phases_s = np.array([kk_phase(np.deg2rad(phi), angles, reflectivities_s) for phi in angles])

delta = phases_p-phases_s

# {{{ Plotten

fig = plt.figure()
ax = fig.add_subplot(211)
plt.plot(angles, reflectivities_p, 'r-')
plt.plot(angles, reflectivities_s, 'b-')
ax = fig.add_subplot(212)
plt.plot(angles, 4*np.rad2deg(delta), 'g-')
fig.savefig('phase.pdf')

# }}}

# vim: foldmethod=marker
