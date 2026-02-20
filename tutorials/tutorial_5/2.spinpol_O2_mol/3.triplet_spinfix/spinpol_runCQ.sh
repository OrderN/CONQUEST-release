#!/bin/bash

# Set env variables
source ../../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Just run Conquest
$CQCMD

# Extract data from Conquest output
#./optpos_extCQ.sh

# Plot data using a python script
#python optpos_pltCQ.py

# Clean around but keep Conquest_out and trajectory.xsf
#./optpos_clnCQ.sh

# Launch vmd using a visualisation data file 
#$VMD_CMD trajectory.xsf -e optpos_visCQ.vmd
