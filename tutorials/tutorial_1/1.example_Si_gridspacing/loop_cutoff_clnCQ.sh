#!/bin/bash

# List of cutoff energy values
values=(50 75 100 120 150)

rm_files=("Conquest_out_ase" 
"Conquest_warnings"
"conquest.bib"
"coord_next.dat"
"eigenvalues.dat"
"hilbert_make_blk.dat"
"InfoGlobal.i00"
"input.log"
"Kmatrix*")

for val in "${values[@]}"; do
    folder=$(printf "gs_%03d" "$val")
    if [ -d $folder ]; then
       cd $folder
       for f in "${rm_files[@]}"; do
          rm -f -- $f 
       done
       cd .. 
    fi
done
