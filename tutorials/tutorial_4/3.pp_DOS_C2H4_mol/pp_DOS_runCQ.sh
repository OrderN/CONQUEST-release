#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Setup PostProcessing command
CQPPCMD="$CQROOT/PostProcessCQ"
echo "PostProcessing cmd: $CQPPCMD"

# Just run Conquest
$CQCMD > Conquest_out_dos

# Just run PostProcessing 
$CQPPCMD > Conquest_out_pp

# Launch xmgrace
#xmgrace -nxy DOS.dat

