# -*- coding: utf-8 -*-

import numpy as np
import polarizer as pol
import matplotlib.pyplot as plt

# init layers
vac = pol.Layer( 'vacuum' )
B4C = pol.Layer( 'B4C', np.complex( 1-0.0850091055, -0.018499583 ), 3e-9 )
Mo = pol.Layer( 'Mo', np.complex( 1-0.2158086, -0.12803553 ), 1e-7 )
sub = pol.Layer( 'SiO2', np.complex( 1-0.0544616096, -0.0341861285 ) )

# build mirror
mirror = pol.Structure( [vac, B4C, Mo, sub] )
mirror.show()

# calculate and plot
angles = np.arange(0, np.pi/2, .01)
r_s = [ mirror.fresnel( angle = a, energy = 60, polarisation = 's' ) for a in angles ]
R_s = [ abs(r)**2 for r in r_s ]

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot( angles, R_s )
fig.savefig('test.pdf')
