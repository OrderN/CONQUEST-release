#!/bin/bash

# Set env variables
source ../../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Just run Conquest
$CQCMD

# Extract data from Conquest output
./spinpol_extCQ.sh
