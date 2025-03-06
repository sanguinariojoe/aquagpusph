#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* fluid.* log.* set0.*.dat set1.*.dat rmse.dat

# Run the "simulation", and clear everything but the main result
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i main.xml

# Check the result
python @TEST_DEST_DIR@/check.py
