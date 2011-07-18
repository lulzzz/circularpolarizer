import polarizer
import numpy as np
import matplotlib.pyplot as plt

# SUBROUTINE polarizer(idx, thickn, n, angles, nA, energy,&
# R_tot, dop, phi)


idx = [np.complex(1, 0), np.complex(.915, -.0185),
        np.complex(.784, -.128), np.complex(.9455, -.0342)]
thickness = [1, 3e-9, 50e-9, 1]
angles = np.arange(np.pi / 4, np.pi / 2, .01)
energy = 60
R, DoP, phi = polarizer.polarizer(idx, thickness, angles, energy)

fig = plt.figure()
ax1 = fig.add_subplot(111)

angles = np.rad2deg(angles)

ax1.plot(angles, DoP)
ax1.plot(angles, R)

# ax2 = ax1.twinx()
# ax2.plot(angles, np.rad2deg(phi), color='red')

# fig.savefig('test.pdf')
fig.show()
