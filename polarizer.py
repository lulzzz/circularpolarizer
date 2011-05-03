# -*- coding: utf-8 -*-

import numpy as np

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
# calculate the imaginary part of reflectivity at angle phi
# kk_real is expected to be [[angles], [reflectivities]]
kk_imag = lambda phi, kk_real: -2/np.pi*sum([phi*kk_real[1]/(kk_real[0]**2-phi**2)])

# }}}


# vim: foldmethod=marker
