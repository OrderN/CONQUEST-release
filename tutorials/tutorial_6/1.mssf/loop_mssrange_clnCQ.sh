#!/bin/bash

# List of cutoff energy values
values=(1.0 5.0 8.0 11.0 14.0)

rm_files=("Conquest_out_ase" 
"Conquest_warnings"
"conquest.bib"
"coord_next.dat"
"eigenvalues.dat"
"hilbert_make_blk.dat"
"InfoGlobal.i00"
"input.log"
"Kmatrix*"
"SFcoeffmatrix*"
)
for val in "${values[@]}"; do
    folder=$(printf "msr_%04.1f" "$val")
    if [ -d $folder ]; then
       cd $folder
       for f in "${rm_files[@]}"; do
          rm -f -- $f 
       done
       cd .. 
    fi
done
