#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* out*.dat

# Run the "simulation" on serial and parallel
/home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -d 2 -i main_serial.xml
mpirun --hostfile hostfile /home/yanez/src/c++/branch_yanez/build/tests/AQUAgpusph -l 00 -d 2 -i main_mpi.xml

# Check the result
python /home/yanez/src/c++/branch_yanez/build/tests/2D/MPI_plane/check.py
