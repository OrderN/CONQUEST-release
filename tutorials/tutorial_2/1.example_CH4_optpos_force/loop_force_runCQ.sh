#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of force tolerance values to test (in Ha/Bohr)
values=(0.0001 0.005 0.001 0.01)

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each force tolerance value
for f in "${values[@]}"; do
    folder="$(printf "f_%.4f" "$f")"
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
    cp C.ion H.ion Conquest_input coords.dat $folder
    cd $folder

    # Generate Conquest_input file with current k-points value
    echo "atommove.maxforcetol ${f}"  >> Conquest_input

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
