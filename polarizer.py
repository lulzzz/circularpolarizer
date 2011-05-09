# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import e,c,h

ii = np.complex(0, 1)

# {{{ single layer-class
class Layer(object):
    def __init__(self, name, index = np.complex(1,0), thickness = -1 ):
        self.name = name
        self.index = index
        self.thickness = thickness

    def _vis(self):
        tw = 45             # textwidth
        text = '%s: N=%1.2f%+1.2fi d=%dnm' % ( self.name, self.index.real, self.index.imag, self.thickness*1e9 )
        return '  |    ' + text + (tw - len(text) - 4) * ' ' + '|\n  |' + tw * '_' + '|\n'
# }}}

# {{{ class for whole structure
class Structure(object):
    def __init__( self, stack ):
        self.stack = stack

    # {{{ generate simple visualization of compound structure
    def show(self):
        # generate simple structure visualization
        text = '\n\n    ambient: %s\n    N=%1.2f%+1.2fi\n' \
                % ( self.stack[0].name, self.stack[0].index.real, self.stack[0].index.imag )
        text += '   ' + 45 * '_' + '\n'
        for layer in self.stack[1:-1]:
            text += layer._vis()
        text += '  |    substrate: %s' % self.stack[-1].name + \
                (30-len(self.stack[-1].name) ) * ' ' + \
                '|\n  |    N=%1.2f %+1.2fi' % ( self.stack[-1].index.real, self.stack[-1].index.imag ) +\
                (len( str(self.stack[-1].index) ) ) * ' ' + '|\n' +\
                '  |' + 45 * ' ' + '|\n\n\n'
        print text
        return text
    # }}}

    # the actual 'worker-function'
    def fresnel( self, angle, energy, polarisation ):
        M = self.__M_ges( self.stack, angle, energy, polarisation )
        return self.__r(M)

    # phase factor from Fresnel equations
    def __b_factor(self, d, n, angle, l):
        return 2 * np.pi * d * n * np.cos(angle) / l

    # calculate refracted angle, Snells law
    def __refracted(self, angle, n1, n2):
        return np.arcsin( np.sin(angle) * n1/ n2 )

    # convert E[eV] --> λ[m]
    def __eV2nm(self, E):
        return h * c / ( e * E )

    # Fresnel-reflectivity from M_ges
    def __r(self, M):
        return M[ 1, 0 ] / M[ 0, 0 ]

    # {{{ matrices
    # matrix for ambient layer, p- and s-polarization
    def __M_amb(self, angle, n, polarisation):
        if polarisation == 'p':
            return .5 * np.matrix( [ [ 1, np.cos(angle)/n ], [ -1, np.cos(angle)/n ] ] )
        elif polarisation == 's':
            return .5 * np.matrix( [ [ 1, 1 / ( np.cos(angle) * n ) ], [ 1, -1 / ( np.cos(angle) * n ) ] ] )

    # any middle layer
    def __M_lay(self, d, n, angle_n, l, polarisation):
        b = self.__b_factor( d, n, angle_n, l )
        if polarisation == 'p':
            M = np.matrix( \
                    [ [ np.cos(b), ii * np.cos(angle_n) * np.sin(b) / n ],\
                    [ ii * n * np.sin(b) / np.cos(angle_n), np.cos(b) ] ] )
            return M
        elif polarisation == 's':
            M = np.matrix( \
                    [ [ np.cos(b), ii * np.sin(b) / ( n * np.cos(angle_n) ) ],\
                    [ ii * n * np.sin(b) * np.cos(angle_n), np.cos(b) ] ] )
            return M

    # matrices for an infinite half-space ( =substrate )
    def __M_sub(self, angle_n, n, polarisation):
        if polarisation == 'p':
            return np.matrix( [ [ np.cos(angle_n) / n, 0 ], [ 1, 0 ] ] )
        elif polarisation == 's':
            return np.matrix( [ [ 1 / ( np.cos(angle_n) * n ), 0 ], [ 1, 0 ] ] )

    # M_ges = M_amb * M_lay1 * ... * M_layN * M_sub
    def __M_ges(self, stack, angle, energy, polarisation):
        ambient, substrate = self.stack[0], self.stack[-1]
        M = self.__M_amb( angle, ambient.index, polarisation )               # starting with vacuum
        l_0 = self.__eV2nm(energy)
        n_old = ambient.index
        for layer in self.stack[1:-1]:                             # going through the finite-thickness layers
            d, n_new = layer.thickness, layer.index
            l = l_0 / n_new
            angle = self.__refracted(angle, n_old, n_new)
            M = M * self.__M_lay(layer.thickness, n_new, angle, l, polarisation)
            n_old = n_new
        n_new = substrate.index                                   # last step: infintie-half-space substrate
        angle = self.__refracted(angle, n_old, n_new)
        M = M * self.__M_sub( angle, n_new, polarisation )
        return M
    # }}}
# }}}

# vim: foldmethod=marker
