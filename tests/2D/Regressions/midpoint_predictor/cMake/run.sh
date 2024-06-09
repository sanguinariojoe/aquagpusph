#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* output.* log.* *pre_predictor* *post_predictor*

# Run the "simulation", and clear everything but the main result
@BINARY_DIR@/AQUAgpusph -d 2 -i main.xml

# Check the result
python @TEST_DEST_DIR@/check.py
