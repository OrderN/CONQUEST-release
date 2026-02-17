#!/bin/bash

# Set env variables
source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# List of basis sets
bases=("SZ" "SZP" "DZP" "TZTP")

# Start timer
start_time="$(date -u +%s.%N)"

# Loop over each basis set
for basis in "${bases[@]}"; do
    folder="bs_$basis"
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi
    # Copy the general input files
    cp coords.dat Conquest_input $folder
    cp "Si_${basis}.ion" $folder
    cd $folder

    # Generate Conquest_input file with current basis
    echo "%block ChemicalSpeciesLabel"    >> Conquest_input
    echo "1   28.085 Si_${basis}"         >> Conquest_input
    echo "%endblock ChemicalSpeciesLabel" >> Conquest_input

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
