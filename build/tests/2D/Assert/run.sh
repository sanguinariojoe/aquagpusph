#!/bin/bash

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat

# Run the "simulation"
/home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -d 2 -i main.xml

set -e

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/2D/Assert/check.py
