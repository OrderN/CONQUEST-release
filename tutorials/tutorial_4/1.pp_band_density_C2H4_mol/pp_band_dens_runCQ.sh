#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Setup PostProcessing command
CQPPCMD="$CQROOT/PostProcessCQ"
echo "PostProc cmd: $CQPPCMD"

# Just run Conquest
$CQCMD

# Just run PostProcessing 
$CQPPCMD > Conquest_pp_out 

# Launch vesta 
$VESTA_CMD *cube 

