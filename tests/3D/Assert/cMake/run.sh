#!/bin/bash

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat

# Run the "simulation"
@BINARY_DIR@/@BINARY_NAME@ -i main.xml

set -e

# Check the result
python @TEST_DEST_DIR@/check.py
