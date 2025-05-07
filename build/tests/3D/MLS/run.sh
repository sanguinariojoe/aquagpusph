#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* fluid.* log.* set0.*.dat set1.*.dat rmse.dat

# Run the "simulation", and clear everything but the main result
/home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -i main.xml

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/3D/MLS/check.py
