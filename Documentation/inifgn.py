#!/usr/bin/env python3

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from numpy import ma

with netCDF4.Dataset("restart.nc") as restart:
    xprimu = restart.variables["xprimu"][:-1]
    xprimv = restart.variables["xprimv"][:-1]

iim = xprimu.size
h = 2 * np.pi / iim

unsddu = 1 / np.sqrt(xprimu)
unsddv = 1 / np.sqrt(xprimv)

deriv_u = np.asmatrix(np.diag(- unsddu * unsddv))
deriv_u[-1, 0] = unsddu[-1] * unsddv[0]

for i in range(iim - 1):
    deriv_u[i, i + 1] = unsddu[i] * unsddv[i + 1]

deriv_v = - deriv_u.T
Delta_u = deriv_u * deriv_v
Delta_v = deriv_v * deriv_u

for Delta in [Delta_u, Delta_v]:
    K = - h**2 * Delta

    # Off-diagonal elements:
    off_diag = ma.masked_array(np.asarray(K), mask = np.identity(iim))

    # Gershgorin radii:
    radii = ma.compressed(np.sum(abs(off_diag), axis = 1))

    plt.matshow(K)
    plt.colorbar()

    plt.figure()
    plt.plot(np.diag(K) - radii)
    plt.ylabel("lower bound of Gershgorin interval")

plt.hot()
plt.show()
