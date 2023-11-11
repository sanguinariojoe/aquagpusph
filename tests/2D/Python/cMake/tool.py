import numpy as np
import os.path as path
import aquagpusph as aqua


def is_approx(v1, v2, tol=1e-5):
    return abs(v1 - v2) < tol


# Read the experimental data
X, Y = np.loadtxt('particles.dat', delimiter=' ', skiprows=0, unpack=True)

def main():
    n = aqua.get("N")
    r = aqua.get("r", 0, n)
    for i in range(n):
        assert(is_approx(r[i][0], X[i]))
        assert(is_approx(r[i][1], Y[i]))

    i = aqua.get("iter")
    i += 1
    aqua.set("iter", i)

    return True
