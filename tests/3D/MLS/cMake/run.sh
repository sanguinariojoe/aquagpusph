#!/bin/bash

set -e

# Run the "simulation", and clear everything but the main result
@BINARY_DIR@/AQUAgpusph -i main.xml
rm -f AQUAgpusph.save.* fluid.* log.*

# Collect the result
rmse=`sed -n '2{p;q;}' rmse.dat`

# We have 500 particles, from which 3 ones have an error of 1, while the other
# ones have errors of ~0. Therefore rmse should be ~sqrt(3 / 500) = 0.077
# Let's get as an acceptable value a 10% bigger error, i.e. 0.085
if (( rmse > 0.085 )); then
    exit 1
fi
