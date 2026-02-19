#!/bin/bash

# Set env variables
source ../../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

export OMP_NUM_THREADS=1


# Start timer
start_time="$(date -u +%s.%N)"

# Just run Conquest
$CQCMD

# End timer
end_time="$(date -u +%s.%N)"

# Print elapsed time
elapsed="$(bc <<<"$end_time-$start_time")"
printf "total of %.3f seconds elapsed\n" "$elapsed"

# Extract data from Conquest output

# Plot data using a python script
python md_analysis.py -s md.stats

# Clean around but keep Conquest_out/trajectory.xsf 
# and md files including velocity.dat
#./md_clnCQ.sh

# Launch vmd using a visualisation data file 
$VMD_CMD trajectory.xsf -e md_visCQ.vmd
