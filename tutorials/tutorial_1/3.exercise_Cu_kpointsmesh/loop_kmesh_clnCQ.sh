#!/bin/bash

# List of kpoints
values=(2 4 6 8)

rm_files=("Conquest_out_ase" 
"Conquest_warnings"
"conquest.bib"
"coord_next.dat"
"eigenvalues.dat"
"hilbert_make_blk.dat"
"InfoGlobal.i00"
"input.log"
"Kmatrix*")

for kp in "${values[@]}"; do
    folder="kp_${kp}x${kp}x${kp}"
    if [ -d $folder ]; then
       cd $folder
       for f in "${rm_files[@]}"; do
          rm -f -- $f 
       done
       cd .. 
    fi
done
