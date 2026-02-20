#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Ion files
ions=("Cu.ion")

# List of temperature values in Ha
values=(0.0001 0.0005 0.001 0.01)

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each k-point triplet
for sigma in "${values[@]}"; do
    folder="$(printf "mp_%.4f" "$sigma")"
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
    cp Cu.ion coords.dat Conquest_input $folder
    cd $folder

    # Generate Conquest_input file with current k-points value
    echo "diag.smearingType 1" >> Conquest_input
    echo "diag.mporder      3" >> Conquest_input
    echo "diag.kt ${sigma}"    >> Conquest_input

    # Run Conquest in the folder
    echo -n "Conquest is running in $folder... "
    $CQCMD > Conquest_out
    echo "done" 
    cd ..
done

# End timer
end_time="$(date -u +%s.%N)"

# Print elapsed time
elapsed="$(bc <<<"$end_time-$start_time")"
printf "total of %.3f seconds elapsed\n" "$elapsed"
