#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Just run Conquest
$CQCMD

# Extract data from Conquest output
./optcell_extCQ.sh

# Plot data using a python script
python optcell_pltCQ.py

# Clean around but keep Conquest_out and trajectory.xsf
./optcell_clnCQ.sh

