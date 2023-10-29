#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat set1.*.dat

# Run the "simulation"
@BINARY_DIR@/AQUAgpusph -i main.xml
