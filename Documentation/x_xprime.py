#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def heavyside(x):
    return np.where(x > 0, 1, 0)

def x(xprime):
    """equation (4) de "Un zoom de deuxième génération pour LMD-Z"."""
    return xprime - 2 * np.pi * heavyside(xprime - np.pi) \
        + 2 * np.pi * heavyside(- xprime - np.pi)

plt.rcParams["axes.labelsize"] = "xx-large"
xprime = np.linspace(-10,10, 1000)
plt.plot(xprime, x(xprime))
plt.ylabel(r"$x$")
plt.xlabel(r"$x'$")
plt.xticks(np.linspace(- 2 * np.pi, 2 * np.pi, 5), 
           (r"$-2 \pi$", r"$- \pi$", "0", r"$\pi$", r"$2 \pi$"))
plt.yticks(np.linspace(- np.pi, np.pi, 3), (r"$- \pi$", "0", r"$\pi$"))

plt.show()
