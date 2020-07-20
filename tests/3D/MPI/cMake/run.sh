#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* out_*.dat

# Run the "simulation"
mpirun --hostfile hostfile @BINARY_DIR@/AQUAgpusph -i main.xml

# Check the result
python @TEST_DEST_DIR@/check.py
