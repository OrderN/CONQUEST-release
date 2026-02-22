#!/bin/bash

# Set env variables
source ../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of files to copy including the block
files=("Si.ion" "coords.dat" "Conquest_template" "block.txt")

# List of kpoints
values=(1 2 4 6)

# Folder basename
base="kp_%03d%03d%03d"

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each cutoff energy value
for val in "${values[@]}"; do
    folder=$(printf "$base" "$val" "$val" "$val")
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy required input files into the folder
cp "${files[@]}" $folder
    cd $folder

    # Generate Conquest_input file with current block
    cp Conquest_template Conquest_input 
    while IFS= read -r line; do
       eval "echo \"$line\"" >> Conquest_input
    done < block.txt

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
