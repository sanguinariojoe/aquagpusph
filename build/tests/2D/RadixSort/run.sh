#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat vars.out

# Run the "simulation"
/home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -d 2 -i main.xml

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/2D/RadixSort/check.py
