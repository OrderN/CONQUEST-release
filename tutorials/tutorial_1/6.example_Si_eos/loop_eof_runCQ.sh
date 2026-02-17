#!/bin/bash

source ../../environment_variables.sh

# Setup Conquest command
CQCMD="$CQMPI $CQNUM $CQROOT/Conquest"
echo "Conquest cmd: $CQCMD"

# Reference lattice constant (in angstroms)
a0_ang=5.4951
ang_to_bohr=1.88973

# Scaling factors: from -5% to +7% around a0
scales=(0.95 0.98 1.00 1.02 1.05 1.07)

# Start timer
start_time="$(date -u +%s.%N)"

for s in "${scales[@]}"; do
    # Compute new lattice constant
    a_bohr=$(echo "$a0_ang * $ang_to_bohr * $s" | bc -l)
    a_ang=$(echo "$a0_ang * $s" | bc -l)

    # Create a folder for this scaling
    folder="$(printf "a_%4.3f" "$a_ang")"
    if [ ! -d $folder ]; then
       mkdir -p $folder
    fi

    # Copy required input files into the folder
    cp Si.ion Conquest_input "$folder/"
    cd $folder

    # Generate the input structure file with updated lattice constant
    cat > "./coords.dat" << EOF
$a_bohr     0.000000    0.000000
0.000000    $a_bohr     0.000000    
0.000000    0.000000    $a_bohr    
8
    0.75000000    0.75000000    0.25000000   1 T T T
    0.00000000    0.50000000    0.50000000   1 T T T
    0.75000000    0.25000000    0.75000000   1 T T T
    0.00000000    0.00000000    0.00000000   1 T T T
    0.25000000    0.75000000    0.75000000   1 T T T
    0.50000000    0.50000000    0.00000000   1 T T T
    0.25000000    0.25000000    0.25000000   1 T T T
    0.50000000    0.00000000    0.50000000   1 T T T
EOF
    # Run Conquest in the folder
    echo -n "Conquest is running in $folder... "
    $CQCMD > Conquest_out
    echo "done" 
    cd ..
done

# End timer
end_time="$(date -u +%s.%N)"

# Compute and print elapsed time
elapsed="$(bc <<<"$end_time-$start_time")"
printf "total of %.3f seconds elapsed\n" "$elapsed"
