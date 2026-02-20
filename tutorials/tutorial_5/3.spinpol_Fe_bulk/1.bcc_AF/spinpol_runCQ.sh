#!/bin/bash

# Set env variables
source ../../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Setup PostProcessing command
CQPPCMD="$CQROOT/PostProcessCQ"
echo "PostProc cmd: $CQPPCMD"

# Just run Conquest
$CQCMD

# Just run Conquest
$CQCMD

# Just run PostProcessing 
$CQPPCMD > Conquest_pp_out

# Extract data from Conquest output
./spinpol_extCQ.sh
