#!/bin/bash

set -e

# Clean up before running, so previous failures are not contaminating the
# results
rm -f AQUAgpusph.save.* log.* set0.*.dat vars.out

# Run the "simulation"
set +e
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i invalid_xml.xml
@BINARY_DIR@/@BINARY_NAME@ -d 2 -i invalid_dom.xml
set -e

# Check the result
python @TEST_DEST_DIR@/check.py
