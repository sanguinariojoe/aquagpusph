#!/usr/bin/env bash

clinfo | grep "Number of devices" > clinfo_platforms
clinfo | grep "Device Type" > clinfo_devices
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python $SCRIPT_DIR/Create.py
chmod +x clean.sh
chmod +x run.sh
./clean.sh
./run.sh
