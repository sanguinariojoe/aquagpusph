# We have 500 particles, from which 3 ones have an error of 1, while the other
# ones have errors of ~0. We add a 5% tolerance
TOL = 1.05 * (3 / 500)**0.5 

with open("rmse.dat", 'r') as f:
    rmse = float(f.readlines()[1])
assert(rmse < TOL)
