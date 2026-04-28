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
$CQCMD

# Just run PostProcessing 
$CQPPCMD > Conquest_pp_out 

grep 'force:   ' Conquest_out | cut -w -f4-6 > tmpfor
tail -8 coords.dat.xsf > tmppos
head -7 coords.dat.xsf > tmphead
paste tmppos tmpfor > tmppos_for
cat tmphead tmppos_for > coords_for.dat.xsf
rm tmppos tmpfor tmphead tmppos_for
