# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e,c,h

#todo: cleanup, classes etc.

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

# {{{ oo stuff comes here

class Layer(object):
    def __init__(self, name, index = np.complex(1,0), thickness = -1 ):
        self.name = name
        self.index = index
        self.thickness = thickness

    def vis( self, textwidth = 35 ):
        text = '%s: N=%1.2f%+1.2fi d=%dnm\n\n' % ( self.name, self.index.real, self.index.imag, self.thickness*1e9 )
        return textwidth * '_' + '\n' + ( textwidth - len(text) ) / 2 * ' ' + text

# class for the whole structure
class Structure(object):
    def __init__( self, stack ):
        self.stack = stack

    def show(self):
        # generate simple structure visualization
        text = '\n\n  Ambient: %s\n  N=%1.2f%+1.2fi\n' % ( self.stack[0].name, self.stack[0].index.real, self.stack[0].index.imag )
        for layer in self.stack[1:-1]:
            text += layer.vis(35)
        text += 35 * '_' + '\n' + \
                '  substrate: %s\n  N=%1.2f %+1.2fi\n' % ( self.stack[-1].name, self.stack[-1].index.real, self.stack[-1].index.imag )
        print text


    def fresnel( self, angle, energy, polarisation ):
        M = M_ges( self.stack, angle, energy, polarisation )
        return r(M)

# {{{ old functions

# closest-match finder
closest = lambda _list, match: min(_list, key = lambda x: abs(x-match) )

nm2eV = lambda l: h * c / ( e * l )
eV2nm = lambda E: h * c / ( e * E )

# phase of complex number
def cplx_phase(c):
    if c.real > 0:
        return np.arctan( c.imag / c.real )
    if c.real < 0 and c.imag >= 0:
        return np.arctan( c.imag / c.real ) + np.pi
    if c.real < 0 and c.imag < 0:
        return np.arctan( c.imag / c.real ) - np.pi

ii = np.complex(0, 1)

b_factor = lambda d, n, angle, l: 2 * np.pi * d * n * np.cos(angle) / l

refracted = lambda phi, n1, n2: np.arcsin( np.sin(phi) * n1/ n2 )

# matrix for ambient layer, p- and s-polarization
# angle of incidence phi is measured from surface-normal!
def M_amb(angle, n, polarisation):
    if polarisation == 'p':
        return .5 * np.matrix( [ [ 1, np.cos(angle)/n ], [ -1, np.cos(angle)/n ] ] )
    elif polarisation == 's':
        return .5 * np.matrix( [ [ 1, 1 / ( np.cos(angle) * n ) ], [ 1, -1 / ( np.cos(angle) * n ) ] ] )
    
# same for arbitrary layer of thickness d, index n, angle of incidence phi (REFRACTED!)
# at wavelength l
def M_lay(d, n, phi_n, l, polarisation):
    b = b_factor( d, n, phi_n, l )
    if polarisation == 'p':
        M = np.matrix( \
                [ [ np.cos(b), ii * np.cos(phi_n) * np.sin(b) / n ],\
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
def M_ges(stack, angle, energy, polarisation):
    ambient = stack[0]
    substrate = stack[-1]
    M = M_amb( angle, ambient.index, polarisation )               # starting with vacuum
    phi_n = angle
    l_0 = eV2nm(energy)
    n_old = ambient.index
    for layer in stack[1:-1]:                             # going through the finite-thickness layers
        d, n_new = layer.thickness, layer.index
        l = l_0 / n_new
        phi_n = refracted(phi_n, n_old, n_new)
        M = M * M_lay(layer.thickness, n_new, phi_n, l, polarisation)
        n_old = n_new
    n_new = substrate.index                                   # last step: infintie-half-space substrate
    phi_n = refracted(phi_n, n_old, n_new)
    M = M * M_sub( phi_n, n_new, polarisation )
    return M

# Fresnel-reflectivity from M_ges
r = lambda M: M[ 1, 0 ] / M[ 0, 0 ]
# }}}

# }}}

# phase_diff = np.array( [ phase_s[i] - p_p for i, p_p in enumerate(phase_p) ] )

# for comparison: data from CXRO-homepage
# deg, cxro_s = np.genfromtxt('/home/dscran/Documents/promotion/circularpolarizer/data/E_60eV_s.dat', unpack=True, skip_header=2, usecols=[0,1])
# deg, cxro_p = np.genfromtxt('/home/dscran/Documents/promotion/circularpolarizer/data/E_60eV_p.dat', unpack=True, skip_header=2, usecols=[0,1])

# fig = plt.figure()
# ax1 = fig.add_subplot( 111 )
# ax1.set_xlabel( u'grazing angle [°]' )
# ax1.set_ylabel( u'reflectivity' )

# plt.plot( np.rad2deg( np.pi/2-angles ), R_s, 'b-', label='s', lw=2 )
# plt.plot( deg, cxro_s, 'bD', label='CXRO' )

# plt.plot( np.rad2deg( np.pi/2-angles ), R_p, 'r-', label='p', lw=2 )
# plt.plot( deg, cxro_p, 'rD', label='CXRO' )

# plt.legend( loc='upper center' )
# 
# ax2 = plt.twinx()
# ax2.set_ylabel( u'phase shift [°]' )

# plt.plot( np.rad2deg( np.pi/2-angles ), -np.rad2deg( phase_s ), 'b.', label='phase_s' )
# plt.plot( np.rad2deg( np.pi/2-angles ), -np.rad2deg( phase_p ), 'r.', label='phase_p' )
# plt.plot( np.rad2deg( np.pi/2-angles ), -np.rad2deg( phase_diff ), 'g-', label='phase shift' )
# plt.legend( loc = 'center right' )

# ax2.set_ylim( ymax=190 )
# ax2.set_yticks( np.arange(0,181,20) )

# plt.legend()

# fig.savefig('/home/dscran/Documents/promotion/circularpolarizer/phase.pdf')

# vim: foldmethod=marker
