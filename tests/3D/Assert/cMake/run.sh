#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat

# Run the "simulation"
if @BINARY_DIR@/AQUAgpusph -i main.xml ; then
    echo "It is supposed that AQUAgpusph is exiting upon an error";
    exit 255;
fi

# Check the result
python @TEST_DEST_DIR@/check.py
