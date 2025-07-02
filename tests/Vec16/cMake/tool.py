import numpy as np
import os.path as path
import aquagpusph as aqua


def is_approx(v1, v2, tol=1e-5):
    return abs(v1 - v2) < tol


# Read the experimental data
data = np.loadtxt('particles.dat', delimiter=' ', skiprows=0, unpack=True)

def main():
    n = aqua.get("N")
    v = aqua.get("v", 0, n)
    for i in range(n):
        for j in range(16):
            assert(is_approx(v[i, j], data[j + 2, i]))

    return True
