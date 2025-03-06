#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* out*.dat

# Run the "simulation" on serial and parallel
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i main_serial.xml
mpirun --hostfile hostfile @BINARY_DIR@/@BINARY_NAME@ -l 00 -d 2 -i main_mpi.xml

# Check the result
python @TEST_DEST_DIR@/check.py
