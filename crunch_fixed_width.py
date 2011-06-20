# -*- coding: utf-8 -*-

import numpy as np
import polarizer as pol
from scipy import interpolate


# {{{ functions
# {{{ imaginary i
ii = np.complex(0, 1)
# }}}


# {{{ closest-match finder
closest = lambda _list, match: min(_list, key=lambda x: abs(x - match))
# }}}


# {{{ get complex refractive index for given material and energy
def get_index(material, energy):
    E, n_r, n_i = np.genfromtxt('data/n_%s.dat' % material,
        skip_header=2, unpack=True)
    sel = [E == closest(E, energy)]
    N = np.complex(1 - n_r[sel], -n_i[sel])
    return N, E[sel][0]
# }}}


# {{{ phase of complex number
def cplx_phase(c):
    if c.real > 0:
        return np.arctan(c.imag / c.real)
    if c.real < 0 and c.imag >= 0:
        return np.arctan(c.imag / c.real) + np.pi
    if c.real < 0 and c.imag < 0:
        return np.arctan(c.imag / c.real) - np.pi
# }}}


# {{{ reflectivity & phasediff from complex fresnel reflection coeff
def fresnel2Rd(r_s, r_p):
    R_s, R_p = abs(r_s) ** 2, abs(r_p) ** 2
    p_s, p_p = [cplx_phase(c) for c in r_s], [cplx_phase(c) for c in r_p]
    delta_p = np.array(p_p) - np.array(p_s)
    return R_s, R_p, delta_p
# }}}


# {{{ return reflectivity for given phase change
def match_refl_to_phase(r_s, r_p, angles, target=np.deg2rad(22.5)):
    R_s, R_p, delta_p = fresnel2Rd(r_s, r_p)
    sel = [delta_p == closest(delta_p, target)]
    return R_s[sel], R_p[sel], angles[sel]
# }}}


# {{{ interpolate calculated phase-shift and reflectivity instead of closest-match
def match_interpolated(r_s, r_p, angles, target=np.deg2rad(22.5)):
    R_s, R_p, delta_p = fresnel2Rd(r_s, r_p)
    # cubic cpline interpolation for phasediff, reflectivities
    delta_p_fit = interpolate.UnivariateSpline(angles, delta_p - target, s=0)
    R_s_fit = interpolate.UnivariateSpline(angles, R_s, s=0)
    R_p_fit = interpolate.UnivariateSpline(angles, R_p, s=0)
    # use Nullstellen (after lowering target angle to zero)
    # there are several Nullstellen; which one to use?
    # don't know, seems to work like this:
    interpolated_angle = delta_p_fit.roots()
    try:
        interpolated_angle = interpolated_angle[1]
    except IndexError:
        interpolated_angle = interpolated_angle[0]
    print np.rad2deg(interpolated_angle)
    interpolated_R_s = R_s_fit(interpolated_angle)
    interpolated_R_p = R_p_fit(interpolated_angle)
    return interpolated_R_s, interpolated_R_p, interpolated_angle
# }}}


# {{{ generate mirror from lists of [material], [thickness]
def gen_mirror(material, thickness, energy):
    mirror = [pol.Layer('vacuum')]
    for i, mat in enumerate(material):
        refr, E = get_index(mat, energy)
        E_diff = abs(E - energy)
        if E_diff > .2:
            print '%s @ %.1feV: target and actual energy differ by %.2feV!'\
                    % (mat, energy, E_diff)
        mirror.append(pol.Layer(mat, refr, 1e-9 * thickness[i]))
    mirror.append(pol.Layer('SiO2', get_index('SiO2', energy)[0]))
    return pol.Structure(mirror)
# }}}


# {{{ cat layer names
def makename(namelist):
    fname = ''
    for s in names:
        fname += s + '-'
    return fname.rstrip('-')
# }}}
# }}}


# tried: Au, B4C/Mo, Si3N4/Mo, BN/Mo, Ru, TiO/Ti,
mnames = [
        ['B4C', 'Mo'],
        # ['Mo'],
        # ['Au'],
        # ['TiO', 'Ti']
        ]
thickness = [
        [3, 300],
        # [300],
        # [300],
        # [5, 300]
        ]
angles = np.arange(0, np.pi, .01)
for nr, names in enumerate(mnames):
    # r = ['# mirror: %s; fmt: energy(eV), R_s(phasedelta=22.5deg), R_p(...)']
    r = []
    fname = makename(names)
    for energy in np.arange(45, 65.1, .3):
        mirror = gen_mirror(names, thickness[nr], energy)
        r_s = np.array([mirror.fresnel(a, energy, 's') for a in angles])
        r_p = np.array([mirror.fresnel(a, energy, 'p') for a in angles])
        # R_s, R_p, angle = match_refl_to_phase(r_s, r_p, angles)
        R_s, R_p, angle = match_interpolated(r_s, r_p, angles)
        print R_s, R_p, angle
        r.append([energy, R_s, R_p, angle])
    np.savetxt('data/interpolated_E-sweep_' + fname + '.txt', r)

# vim: foldmethod=marker
