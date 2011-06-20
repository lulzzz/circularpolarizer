# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e, c, h


# {{{ single layer-class
class Layer(object):
    def __init__(self, name, index=np.complex(1, 0), d=-1):
        self.name = name
        self.index = index
        self.d = d

    def _vis(self):
        text_d = '' if self.d == -1 else ', d=%dnm' % (self.d * 1e9)
        text = '%s: N=%1.2f%+1.2fi%s' \
                % (self.name, self.index.real, self.index.imag, text_d)
        return text
# }}}


# {{{ class for whole structure
class Structure(object):
    # convenient, isn't it?
    ii = np.complex(0, 1)

    def __init__(self, stack):
        self.stack = stack

    # {{{ generate simple visualization of compound structure
    def show(self):
        text = ''
        for layer in self.stack:
            text += layer._vis() + ' | '
        return text
    # }}}

    # the actual 'worker-function'
    def fresnel(self, angle, energy, polarisation):
        M = self.__M_ges(self.stack, angle, energy, polarisation)
        return self.__r(M)

    # {{{ helpers for reflectivity calculation
    # phase factor from Fresnel equations
    def __b_factor(self, d, n, angle, l):
        return 2 * np.pi * d * n * np.cos(angle) / l

    # calculate refracted angle, Snells law
    def __refracted(self, angle, n1, n2):
        return np.arcsin(np.sin(angle) * n1 / n2)

    # convert E[eV] --> λ[m]
    def __eV2nm(self, E):
        return h * c / (e * E)

    # Fresnel-reflectivity from M_ges
    def __r(self, M):
        return M[1, 0] / M[0, 0]

    # {{{ matrices
    # matrix for ambient layer, p- and s-polarization
    def __M_amb(self, angle, n, polarisation):
        if polarisation == 'p':
            return .5 * np.matrix([[1, np.cos(angle) / n],
                [-1, np.cos(angle) / n]])
        elif polarisation == 's':
            return .5 * np.matrix([[1, 1 / (np.cos(angle) * n)],
                [1, -1 / (np.cos(angle) * n)]])

    # any middle layer
    def __M_lay(self, d, n, angle_n, l, polarisation):
        b = self.__b_factor(d, n, angle_n, l)
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
    def __M_sub(self, angle_n, n, polarisation):
        if polarisation == 'p':
            return np.matrix([[np.cos(angle_n) / n, 0], [1, 0]])
        elif polarisation == 's':
            return np.matrix([[1 / (np.cos(angle_n) * n), 0], [1, 0]])

    # M_ges = M_amb * M_lay1 * ... * M_layN * M_sub
    def __M_ges(self, stack, angle, energy, polarisation):
        ambient, substrate = self.stack[0], self.stack[-1]
        # starting with vacuum
        M = self.__M_amb(angle, ambient.index, polarisation)
        l_0 = self.__eV2nm(energy)
        n_old = ambient.index
        # going through the finite-thickness layers
        for layer in self.stack[1:-1]:
            d, n_new = layer.d, layer.index
            l = l_0 / n_new
            angle = self.__refracted(angle, n_old, n_new)
            M = M * self.__M_lay(layer.d, n_new,
                    angle, l, polarisation)
            n_old = n_new
        # last step: infintie-half-space substrate
        n_new = substrate.index
        angle = self.__refracted(angle, n_old, n_new)
        M = M * self.__M_sub(angle, n_new, polarisation)
        return M
    # }}}
    # }}}


'''
this is a class for simple geometric visualization, useful to show the path of
the lightbeam along the 4 mirrors of the polarizer.  A straightline-object is
provided, that stores lines in basepoint-direction-form. Direction is given as
angle.
Methods are 
> straightline.isct(otherline)
that returns the intersection point of 'straightline' and 'otherline'.
'otherline' has to be a straighline object too.
> straightline.mirror(mirrorline) returns intersection point
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
        raypath.append(ray.goto(mirrorlist[-1].bp[0] + 50 - bp[0]))
        return np.array(raypath)

# }}}

# vim: foldmethod=marker
