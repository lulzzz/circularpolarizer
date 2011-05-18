# -*- coding: utf-8 -*-

import numpy as np
import polarizer as pol
import matplotlib.pyplot as plt
from scipy.constants import h, e, c
from itertools import product

# {{{ helper functions

def cartesian(arrays, out=None):
    """Generate a cartesian product of input arrays."""

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

# import material data
def import_constants(filename):
    n_real, n_imag = np.genfromtxt( filename, skip_header = 3, unpack = True )
    return n_real, n_imag

# imaginary i
ii = np.complex(0, 1)

# closest-match finder
closest = lambda _list, match: min(_list, key = lambda x: abs(x-match) )

# get complex refractive index for given _material_ at _energy_
def get_index( material, energy):
    E, n_r, n_i = np.genfromtxt( 'data/n_%s.dat' % material, skip_header = 2, unpack = True )
    sel = [ E == closest(E, energy) ]
    E = E[sel]
    N = np.complex( 1-n_r[sel], -n_i[sel] )
    return E, N

# phase of complex number
def cplx_phase(c):
    if c.real > 0:
        return np.arctan( c.imag / c.real )
    if c.real < 0 and c.imag >= 0:
        return np.arctan( c.imag / c.real ) + np.pi
    if c.real < 0 and c.imag < 0:
        return np.arctan( c.imag / c.real ) - np.pi

def match_refl( list_refl, list_phase, target = 22.5 ):
    return list_refl[ list_phase == closest( list_phase, 22.5 ) ]
# }}}

# i'd like to play through some material/thickness combinations:
# single-layer müssten extra gehandhabt werden, also einfach zwei gleichartige nehmen...
crunch = ( \
	( ( 'B4C'   , np.arange( 1 , 10, 3 ) )   , ( 'Mo', np.arange( 10, 51, 5 ) ) ),\
	( ( 'Al2O3' , np.arange( 1 , 10, 3 ) )   , ( 'Al', np.arange( 10, 51, 5 ) ) ),\
	( ( 'TiO'   , np.arange( 1 , 10, 3 ) )   , ( 'Ti', np.arange( 10, 51, 5 ) ) ),\
	( ( 'Au'    , np.arange( 1, 5, 5 ) )   , ( 'Au', np.arange( 10, 51, 5 ) ) ),\
	( ( 'Pt'    , np.arange( 1, 5, 5 ) )   , ( 'Pt', np.arange( 10, 51, 5 ) ) ),\
	( ( 'SiC'   , np.arange( 1, 5, 5 ) )   , ( 'SiC', np.arange( 10, 51, 5 ) ) ),\
	)

# energy in eV
energy = 60
# scanning angles
angles = np.arange(0, np.pi/2, .01)

fig = plt.figure(figsize = (16,6) )

# initialize ambient and substrate 'layers'
vac = pol.Layer( 'vacuum' )
E, N = get_index( 'SiO2', energy )
sub = pol.Layer( 'SiO2', N )

for nr, mat_pair in enumerate(crunch):
    list_Rp = []
    list_Rs = []
    nicename = '' 
    grid = []
    namelist = []
    index_of_refraction = []
    ax = fig.add_subplot( 1, len(crunch), nr+1 )
    for layer in mat_pair:
	name = layer[0]
	grid = layer[1] if not any(grid) else cartesian( [grid, layer[1]] )
	namelist.append( name )
	nicename += name + '/'
	index_of_refraction.append( get_index( name, energy )[1] )
    for comb in grid:
	mirror = [ vac ]
	try:
	    for i, thickness in enumerate(comb):
		mirror.append( pol.Layer( namelist[i], index_of_refraction[i], thickness*1e-9 ) )
	except TypeError:
	    # we're probably having only one layer, so iteration fails
	    mirror.append( pol.Layer( namelist[0], index_of_refraction[0], comb*1e-9 ) )
	mirror.append( sub )
	mirror = pol.Structure( mirror )

        r_s = np.array( [ mirror.fresnel( angle = a, energy = energy, polarisation = 's' ) for a in angles ] )
        R_s = np.array( [ abs(r)**2 for r in r_s ] )
        phase_s = np.array( [ cplx_phase(r) for r in r_s ] )

        r_p = np.array( [ mirror.fresnel( angle = a, energy = 60, polarisation = 'p' ) for a in angles ] )
        R_p = np.array( [ abs(r)**2 for r in r_p ] )
        phase_p = np.array( [ cplx_phase(r) for r in r_p ] )
        phase_delta = phase_p - phase_s

	# this is the important value: absolute reflectivity for phase-shift of 90°/4
	R_p = R_p[ phase_delta == closest( phase_delta, np.deg2rad(22.5) ) ] 
	R_s = R_s[ phase_delta == closest( phase_delta, np.deg2rad(22.5) ) ] 
	list_Rp.append([ comb, R_p ])
	list_Rs.append([ comb, R_s ])

	try:
	    pos = comb[1]
	    pltcol = lambda val, maxval: float(maxval)/val
	    m, v = float(max(grid[:,0])), float(comb[0])
	    c = ( (m-v)/(m+v), 1-(m-v)/(m+v), 0 )
	except:
	    pos = comb
	    pltcol = (0, 0, 1)
	plt.plot( pos, R_s**4, 'o', color=c )

    plt.title( nicename.rstrip('/') )
    ax.set_ylim( 0, .20 )
    ax.set_xticks( np.arange(10,51,10) )

fig.subplots_adjust( left=.05, right=.98, wspace=.25 )
fig.savefig('/home/dscran/documents/promotion/circularpolarizer/test_crunch.pdf')

# vim: foldmethod=marker
