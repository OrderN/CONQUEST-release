#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Setup PostProcessing command
CQPPCMD="$CQROOT/PostProcessCQ"
echo "PostProc cmd: $CQPPCMD"

# Run Conquest for SCF
cp Conquest_input_scf Conquest_input
$CQCMD > Conquest_out_scf

# Run Conquest for non-SCF and bands 
cp Conquest_input_bnd Conquest_input
$CQCMD > Conquest_out_bnd

# Just run PostProcessing 
$CQPPCMD > Conquest_pp_out 

