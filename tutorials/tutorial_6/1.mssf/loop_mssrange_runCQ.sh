#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of cutoff energy values in Ha
values=(1.0 5.0 8.0 11.0 14.0)

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each cutoff energy value
for val in "${values[@]}"; do
    folder=$(printf "msr_%04.1f" "$val")
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
    cp Si.ion coords.dat Conquest_input $folder
    cd $folder

    # Generate Conquest_input file with current cutoff value
    echo "%block Si" >> Conquest_input
    echo "Atom.NumberOfSupports 4"  >> Conquest_input
    echo "Atom.MultisiteRange $val" >> Conquest_input
    echo "Atom.LFDRange $val" >> Conquest_input
    echo "%endblock" >> Conquest_input

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
