#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of kpoints
values=(1 2 4 6)

# Loop over each k-points value
for kp in "${values[@]}"; do
    folder="kp_${kp}x${kp}x${kp}"
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
    cp Si.ion coords.dat Conquest_input $folder
    cd $folder

    # Generate Conquest_input file with current k-points value
    echo "diag.mpmesh  True"  >> Conquest_input
    echo "diag.mpmeshx ${kp}" >> Conquest_input
    echo "diag.mpmeshy ${kp}" >> Conquest_input
    echo "diag.mpmeshz ${kp}" >> Conquest_input

    # Run Conquest in the folder
    echo -n "Conquest is running in $folder... "
    $CQCMD > Conquest_out
    echo "done" 
    cd ..
done

