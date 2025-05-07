#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* out_*.dat

# Run the "simulation"
mpirun --hostfile hostfile /home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -l 00 -d 2 -i main.xml

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/2D/MPI_sync/check.py
