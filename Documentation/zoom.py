#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def F(xtild):
    if xtild == 0:
        result = 1
    elif xtild == np.pi:
        result = - 1
    else:
        result = np.tanh(tau * (delta / 2 - xtild) / xtild / (np.pi - xtild))
    return result

try:
    tau = float(input("tau = (default 3) ?"))
except ValueError:
    tau = 3

try:
    delta = float(input("delta = (default 2 pi / 5) ?"))
except ValueError:
    delta = 2 * np.pi / 5

nmax = 30000
xtild = np.linspace(0, np.pi, nmax + 1)
fhyp = np.array([F(xtild_i) for xtild_i in xtild])

try:
    grossismx = float(input("grossismx = (default 2) ?"))
except ValueError:
    grossismx = 2

ffdx = np.empty(nmax + 1)
ffdx[0] = 0
for i in range(nmax):
    F_dxtild, abserr = integrate.quad(F, xtild[i], xtild[i + 1])
    ffdx[i + 1] = ffdx[i] + F_dxtild

beta = (grossismx * ffdx[nmax] - np.pi) / (ffdx[nmax] - np.pi)
print("beta =", beta)

# Plots:

plt.rcParams["axes.labelsize"] = "xx-large"

plt.figure()
plt.xlim(0, np.pi)
plt.plot(xtild, fhyp)
plt.ylabel(r"$F$")
plt.xlabel(r"$\tilde x$")
plt.axvline(delta / 2, linestyle = "--", color = "black")
plt.text(delta / 2, 1.02, r"$\delta / 2$", horizontalalignment = 'center')
plt.axhline(linestyle = "--", color = "black")

plt.figure()
plt.plot(xtild,  beta + (grossismx - beta) * fhyp)
plt.ylabel(r"$g$")
plt.xlabel(r"$\tilde x$")

plt.figure()
plt.axes(aspect = "equal")
plt.xlim(0, np.pi)
plt.ylim(0, np.pi)
plt.plot(xtild,  beta * xtild + (grossismx - beta) * ffdx)
plt.plot(xtild,  xtild, linestyle = "--", color = "black")
plt.ylabel(r"$x_f$")
plt.xlabel(r"$\tilde x$")

plt.show()
