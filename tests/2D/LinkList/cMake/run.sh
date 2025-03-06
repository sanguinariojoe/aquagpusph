#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat vars.out

# Run the "simulation"
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i main.xml

# Check the result
python @TEST_DEST_DIR@/check.py
