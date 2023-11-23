#!/bin/bash

set -e

# Run the "simulation", and clear everything but the main result
@BINARY_DIR@/AQUAgpusph -i main.xml
rm -f AQUAgpusph.save.* fluid.* log.*

# Check the result
python @TEST_DEST_DIR@/check.py
