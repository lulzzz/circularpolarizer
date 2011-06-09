# -*- coding: utf-8 -*-

import numpy as np
import polarizer as pol


# imaginary i
ii = np.complex(0, 1)

# closest-match finder
closest = lambda _list, match: min(_list, key=lambda x: abs(x - match))


# get complex refractive index for given _material_ at _energy_
def get_index(material, energy):
    E, n_r, n_i = np.genfromtxt('data/n_%s.dat' % material,
        skip_header=2, unpack=True)
    sel = [E == closest(E, energy)]
    N = np.complex(1 - n_r[sel], -n_i[sel])
    return N


# phase of complex number
def cplx_phase(c):
    if c.real > 0:
        return np.arctan(c.imag / c.real)
    if c.real < 0 and c.imag >= 0:
        return np.arctan(c.imag / c.real) + np.pi
    if c.real < 0 and c.imag < 0:
        return np.arctan(c.imag / c.real) - np.pi


# return reflectivity for given phase change
def match_refl_to_phase(r_s, r_p, target=np.deg2rad(22.5)):
    R_s, R_p = abs(r_s) ** 2, abs(r_p) ** 2
    p_s, p_p = [cplx_phase(c) for c in r_s], [cplx_phase(c) for c in r_p]
    delta_p = np.array(p_p) - np.array(p_s)
    sel = [delta_p == closest(delta_p, target)]
    R_s, R_p = R_s[sel], R_p[sel]
    return R_s


# generate all possible combinations for a number of input arrays
def combinations(lists):
    r = [[]]
    for x in lists:
        r = [i + [y] for y in x for i in r]
    return r


# input a layer structure
def crunch(names, thickness, parameter):
    energy, angles = parameter
    refr = [get_index(name, energy) for name in names]
    sub = pol.Layer('SiO2', get_index('SiO2', energy))
    r = []
    for comb in combinations(thickness):
        mirror = [pol.Layer('vacuum')]
        for l in [[name, refr[i], comb[i]] for i, name in enumerate(names)]:
            mirror.append(pol.Layer(l[0], l[1], 1e-9 * l[2]))
        mirror.append(sub)
        mirror = pol.Structure(mirror)
        print mirror.show()
        r_s = np.array([mirror.fresnel(a, energy, 's') for a in angles])
        r_p = np.array([mirror.fresnel(a, energy, 'p') for a in angles])
        r += [comb + [match_refl_to_phase(r_s, r_p)]]
    return r


mnames = [['B4C', 'Mo'], ['Rh'], ['Ru'], ['TiO', 'Ti']]
# mnames = [['W'], ['Pd']]
thickness = [
        [np.arange(2, 10, 40), np.arange(10, 100, 100)],
        [np.arange(10, 100, 100)],
        [np.arange(10, 100, 100)],
        [np.arange(2, 10, 40), np.arange(10, 100, 100)]]
for energy in np.arange(55, 65, .1):
    angles = np.arange(0, np.pi, .01)
    for nr, names in enumerate(mnames):
        r = crunch(names, thickness[nr], [energy, angles])
        fname = ''
        for name in names:
            fname += name + '-'
        np.savetxt('data/' + fname.rstrip('-') + '_E=%1.1feV.txt' % energy, r)
