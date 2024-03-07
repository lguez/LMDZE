#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import glob

plt.ylim(0,3)
plt.xticks([])
plt.ylabel("index")

x = 0.

for half_level_file in glob.glob("half_level_*.csv"):
    full_level_file = "full_level_" + half_level_file[11:-4] + ".csv"

    p_half = np.loadtxt(half_level_file, skiprows=1, usecols=(3,))
    p_half = np.append(p_half, 0)

    p_full = np.loadtxt(full_level_file, skiprows=1, usecols=(0,))

    x = x + 0.3

    for i, v in enumerate(p_half):
        plt.text(x, i + .5, str(v))

    for i, v in enumerate(p_full):
        plt.text(x, i + 1, str(v), color="red")

##plt.show()
plt.savefig("test_disvert.pdf")
