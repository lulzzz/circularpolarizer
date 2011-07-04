# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e, c, h


# {{{ class for whole structure
class mirror(object):
    # convenient, isn't it?
    ii = np.complex(0, 1)

    def __init__(self, mirrordef):
        self.names = mirrordef['names']
        self.thickness = 1e-9 * np.array(mirrordef['thickness'])
        self.E = mirrordef['energy']
        try:
            self.index = mirrordef['index']
        except KeyError:
            print 'No refractive indices given, looking for datafiles.'
            self.index = [get_index(n, self.E) for n in self.names]
        try:
            self.ambient = mirrordef['ambient']
        except KeyError:
            self.ambient = np.complex(1,0)
            # set ambient to vacuum if undefined

    # the actual 'worker-function'
    def fresnel(self, angle, polarisation):
        M = self.__M_ges(angle, polarisation)
        return r(M)

    def get_Rd(self, angle):
        r_s = self.fresnel(angle, 's')
        r_p = self.fresnel(angle, 'p')
        return fresnel2Rd(r_s, r_p)


    # {{{ helpers for reflectivity calculation
    # {{{ matrices
    # matrix for ambient layer, p- and s-polarization
    def __M_amb(self, angle, polarisation):
        if polarisation == 'p':
            return .5 * np.matrix([[1, np.cos(angle) / self.ambient],
                [-1, np.cos(angle) / self.ambient]])
        elif polarisation == 's':
            return .5 * np.matrix([[1, 1 / (np.cos(angle) * self.ambient)],
                [1, -1 / (np.cos(angle) * self.ambient)]])

    # any middle layer
    def __M_lay(self, d, n, angle_n, l, polarisation):
        b = b_factor(d, n, angle_n, l)
        if polarisation == 'p':
            M = np.matrix( \
                    [[np.cos(b), self.ii * np.cos(angle_n) * np.sin(b) / n],\
                    [self.ii * n * np.sin(b) / np.cos(angle_n), np.cos(b)]])
            return M
        elif polarisation == 's':
            M = np.matrix(\
                    [[np.cos(b), self.ii * np.sin(b) / (n * np.cos(angle_n))],\
                    [self.ii * n * np.sin(b) * np.cos(angle_n), np.cos(b)]])
            return M

    # matrices for an infinite half-space ( =substrate )
    def __M_sub(self, angle_n, polarisation):
        if polarisation == 'p':
            return np.matrix([[np.cos(angle_n) / self.index[-1], 0], [1, 0]])
        elif polarisation == 's':
            return np.matrix([[1 / (np.cos(angle_n) * self.index[-1]), 0], [1, 0]])

    # M_ges = M_amb * M_lay1 * ... * M_layN * M_sub
    def __M_ges(self, angle, polarisation):
        # starting with vacuum
        M = self.__M_amb(angle, polarisation)
        l_0 = eV2nm(self.E)
        n_old = self.ambient
        index = self.index
        # going through the finite-thickness layers
        for i, layer in enumerate(self.names[:-1]):
            d, n_new = self.thickness[i], index[i]
            l = l_0 / n_new
            angle = refracted(angle, n_old, n_new)
            M = M * self.__M_lay(self.thickness[i], n_new,
                    angle, l, polarisation)
            n_old = n_new
        # last step: infinite-half-space substrate
        n_new = index[-1]
        angle = refracted(angle, n_old, n_new)
        M = M * self.__M_sub(angle, polarisation)
        return M

    def info(self):
        ds = lambda d: '%dnm' % d if d != -1 else 'substrate'
        info = ''.join(['________\n| %s\n| %s\n' % (n, ds(self.thickness[i] * 1e9))
            for i, n in enumerate(self.names)])
        info += '\n\nE=%deV' % self.E
        return info

    # }}}
    # }}}
# }}}

# {{{ straightline class
'''
this is a class for simple geometric visualization, useful to show the path of
the lightbeam along the 4 mirrors of the polarizer.  A straightline-object is
provided, that stores lines in basepoint-direction-form. Direction is given as
angle.
Methods are 
> straightline.isct(otherline)
that returns the intersection point of 'straightline' and 'otherline'.
'otherline' has to be a straighline object too.
> straightline.mirror(mirrorline) 
returns intersection point and angle corresponding to a reflection of
straightline on mirrorline.
Capital letters denote points, small ones stand for direction of straights.
'''
class straightline(object):
    def __init__(self, basepoint, angle):
        self.bp = basepoint
        self.angle = angle
        self.d = np.array([np.cos(angle), np.sin(angle)])

    def goto(self, l):
        return self.bp + l * self.d

    # intersection with other straightline
    def isct(self, line1):
        E0x, E0y = self.bp
        e0x, e0y = self.d
        E1x, E1y = line1.bp
        e1x, e1y = line1.d
        # solution of l for straight==line1
        l = -(e1x * E1y - e1y * E1x + E0x * e1y - E0y * e1x) /\
                (e0x * e1y - e0y * e1x)
        return l

    # reflection on a mirror straight, returns basepoint and angle of
    # reflected straightline
    def refl(self, mirror):
        return self.goto(self.isct(mirror)), 2 * mirror.angle - self.angle

    def mirrorpath(self, mirrorlist):
        ray = self
        raypath = [self.bp]
        for m in mirrorlist:
            bp, angle = ray.refl(m)
            ray = straightline(bp, angle)
            raypath.append(bp)
        raypath.append(ray.goto(mirrorlist[-1].bp[0] + 100 - bp[0]))
        return np.array(raypath)
# }}}

# {{{ reflectivity & phasediff from complex fresnel reflection coeff
def fresnel2Rd(r_s, r_p):
    R_s, R_p = abs(r_s) ** 2, abs(r_p) ** 2
    p_s, p_p = cplx_phase(r_s), cplx_phase(r_p)
    return R_s, R_p, p_p - p_s
# }}}

# {{{ closest match finder
closest = lambda _list, match: min(_list, key=lambda x: abs(x - match))
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

#{{{ get complex refractive index for given material and energy
def get_index(material, energy):
    E, n_r, n_i = np.genfromtxt('data/n_%s.dat' % material,
        skip_header=2, unpack=True)
    sel = [E == closest(E, energy)]
    N = np.complex(1 - n_r[sel], -n_i[sel])
    if energy - E[sel][0] > .5:
        print '%s: target and actual energy differ by %.2deV!'\
                % (material, E[sel][0])
    return N
# }}}

# {{{ calculate refracted angle, Snells law
def refracted(angle, n1, n2):
    return np.arcsin(np.sin(angle) * n1 / n2)
# }}}

# {{{ convert E[eV] --> λ[m]
def eV2nm(E):
    return h * c / (e * E)
# }}}

# {{{ Fresnel-reflectivity from M_ges
def r(M):
    return M[1, 0] / M[0, 0]
# }}}

# {{{ phase factor from Fresnel equations
def b_factor(d, n, angle, l):
    return 2 * np.pi * d * n * np.cos(angle) / l
# }}} 

# vim: foldmethod=marker
