#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return (3 * x**2 - 4 * x + 2) / (2 * x**2 - 2 *x  + 1)

pref = 1013.25
p = np.linspace(0, pref)
plt.plot(f(p / pref), p)
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel("zvert")
plt.ylabel("p (hPa)")
plt.savefig("zvert.pdf")
