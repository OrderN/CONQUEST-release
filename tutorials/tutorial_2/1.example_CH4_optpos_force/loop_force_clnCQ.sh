#!/bin/bash

# List of force tolerance values to test (in Ha/Bohr)
values=(0.0001 0.005 0.001 0.01)

rm_files=("Conquest_out_ase" 
"Conquest_warnings"
"conquest.bib"
"coord_next.dat"
"eigenvalues.dat"
"hilbert_make_blk.dat"
"InfoGlobal.i00"
"input.log"
"Kmatrix*"
"make_prt.dat"
"UpdatedAtoms.dat")

for f in "${values[@]}"; do
    folder=$(printf "f_%.4f" "$f")
    if [ -d $folder ]; then
       cd $folder
       for f in "${rm_files[@]}"; do
          rm -f -- $f 
       done
       cd .. 
    fi
done
