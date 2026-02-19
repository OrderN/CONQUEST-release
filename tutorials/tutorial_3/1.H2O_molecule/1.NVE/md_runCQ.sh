#!/bin/bash

# Set env variables
source ../../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# link NVT directory
ln -s ../0.NVT/reference restart_files

# For a restart we need
cp restart_files/InfoGlobal.i00 .
cp restart_files/Kmatrix* .
cp restart_files/md.checkpoint .
cp restart_files/md.position .

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
#awk 'BEGIN{print "step   pe   ke   H   T   P"}1' md.stats > tmp
#mv tmp md.stats 
#python md_analysis.py -s md.stats

# Clean around but keep Conquest_out/trajectory.xsf and md files
#./md_clnCQ.sh

# Launch vmd using a visualisation data file 
#$VMD_CMD trajectory.xsf -e md_visCQ.vmd
