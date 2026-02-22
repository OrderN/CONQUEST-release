#!/bin/bash

# Set env variables
source ../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of files to copy
files=("Si.ion" "coords.dat" "Conquest_template")

# List of cutoff energy values in Ha
values=(50 75 100 120 150)

# Folder basename
base="gs_%03d"

# String to edit in Conquest_input
string="grid.gridcutoff "

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each cutoff energy value
for val in "${values[@]}"; do
    folder=$(printf "$base" "$val")
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
    cp "${files[@]}" $folder
    cd $folder

    # Generate Conquest_input file with current cutoff value
    cp Conquest_template Conquest_input
    echo $string $val >> Conquest_input 

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
