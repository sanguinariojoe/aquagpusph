#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat

# Run the "simulation"
/home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -i main.xml

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/3D/LinkList/check.py
