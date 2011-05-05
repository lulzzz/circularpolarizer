# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e,c,h

# calculate the complex refraction coefficient (Fresnel) for a multilayer-structure
# matrices are from:
# Harland G. Tompkins, Eugene A. Irene: 'Handbook of ellipsometry' 

# Structure looks like this:
# 
# | Layer 0: Vacuum (n_vac, phi_vac) | Layer 1: e.g. B4C (n1, phi1, b1, d1) | Layer N: (...) | Layer N+1: Substrate (n_sub, phi_sub)
#
# input data should look like:
# ambient = [ phi, n , lambda ]
# layerlist = [ [ d_1, n_1, l_1 ], ... , [ d_N, n_N, l_N ] ]
# substrate = [ n ]

# {{{ definitions

lambda2eV = lambda l: h * c / ( e * l )
eV2lambda = lambda E: h * c / ( e * E )

ii = np.complex(0, 1)

b_factor = lambda d, n, phi, l: 2 * np.pi * d * n * np.cos(phi) / l

refracted = lambda phi, n1, n2: np.arcsin( np.sin(phi) * n1/ n2 )

# matrix for ambient layer, p- and s-polarization
# angle of incidence phi is measured from surface-normal!
def M_amb(phi, n, polarisation):
    if polarisation == 'p':
        return .5 * np.matrix( [ [ 1, np.cos(phi)/n ], [ -1, np.cos(phi)/n ] ] )
    elif polarisation == 's':
        return .5 * np.matrix( [ [ 1, 1 / ( np.cos(phi) * n ) ], [ 1, -1 / ( np.cos(phi) * n ) ] ] )
    
# same for arbitrary layer of thickness d, index n, angle of incidence phi (REFRACTED!)
# at wavelength l
def M_lay(d, n, phi_n, l, polarisation):
    b = b_factor( d, n, phi_n, l )
    if polarisation == 'p':
        M = np.matrix( \
                [ [ np.cos(b), -ii * np.cos(phi_n) * np.sin(b) / n ],\
                [ ii * n * np.sin(b) / np.cos(phi_n), np.cos(b) ] ] )
        return M
    elif polarisation == 's':
        M = np.matrix( \
                [ [ np.cos(b), ii * np.sin(b) / ( n * np.cos(phi_n) ) ],\
                [ ii * n * np.sin(b) * np.cos(phi_n), np.cos(b) ] ] )
        return M

# matrices for an infinite half-space ( =substrate )
def M_sub(phi_n, n, polarisation):
    if polarisation == 'p':
        return np.matrix( [ [ np.cos(phi_n) / n, 0 ], [ 1, 0 ] ] )
    elif polarisation == 's':
        return np.matrix( [ [ 1 / ( np.cos(phi_n) * n ), 0 ], [ 1, 0 ] ] )

# M_ges = M_amb * M_lay1 * ... * M_layN * M_sub
def M_ges(ambient, layerlist, substrate, polarisation):
    M = M_amb( ambient[0], ambient[1], polarisation )               # starting with vacuum
    phi_n, n_old, l_0 = ambient
    for layer in layerlist:                             # going through the finite-thickness layers
        d, n_new = layer
        l = l_0 / n_new
        phi_n = refracted(phi_n, n_old, n_new)
        M = M * M_lay(d, n_new, phi_n, l, polarisation)
        n_old = n_new
    n_new = substrate[0]                                   # last step: infintie-half-space substrate
    phi_n = refracted(phi_n, n_old, n_new)
    M = M * M_sub( phi_n, n_new, polarisation )
    return M

def _test_M_ges(ambient, layerlist, substrate, polarisation):
    M = M_amb( ambient[0], ambient[1], polarisation )               # starting with vacuum
    phi_n, n_old, l_0 = ambient
    for layer in layerlist:                             # going through the finite-thickness layers
        d, n_new = layer
        l = l_0 / n_new
        phi_n = refracted(phi_n, n_old, n_new)
        M = M * M_lay(d, n_new, phi_n, l, polarisation)
        n_old = n_new
    n_new = substrate[0]                                   # last step: infintie-half-space substrate
    phi_n = refracted(phi_n, n_old, n_new)
    M = M * M_sub( phi_n, n_new, polarisation )
    return M

# Fresnel-reflectivity from M_ges
r = lambda M: M[ 1, 0 ] / M[ 0, 0 ]

# phase of complex number
cplx_phase = lambda c: np.arctan( c.imag / c.real )

# }}}

# n_B4C @ 60eV
n_B4C = np.complex( 1-0.0850091055, -0.018499583 )
# n_Mo @ 60eV
n_Mo = np.complex( 1-0.2158086, -0.12803553 )

n_vac = np.complex( 1, 0 )
n_SiO2 = np.complex( 1-0.0544616096, -0.0341861285 )

layerlist = [ [ 3e-9, n_B4C ], [ 100e-9, n_Mo ] ]
substrate = [ n_SiO2 ]
l = eV2lambda(60)

angles = np.arange( 0.01, np.pi/2, .05 )

# shift = [ cplx_phase( r( M_ges( [a, n_vac, l], layerlist, substrate, 's' ) ) ) \
        # - cplx_phase( r( M_ges( [a, n_vac, l], layerlist, substrate, 'p' ) ) ) for a in angles ]

r_s = [ r( M_ges( [a, n_vac, l], layerlist, substrate, 's' ) ).real for a in angles ]
r_p = [ r( M_ges( [a, n_vac, l], layerlist, substrate, 'p' ) ).real for a in angles ] 

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(np.pi/2-angles, r_s, 'r-', label='s')
plt.plot(np.pi/2-angles, r_p, 'g-', label='p')
plt.legend()
fig.savefig('/home/dscran/Documents/promotion/circularpolarizer/phase.pdf')
# vim: foldmethod=marker
