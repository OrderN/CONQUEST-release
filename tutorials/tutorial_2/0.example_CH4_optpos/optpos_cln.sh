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
"Kmatrix*"
"trajectory.xsf"
"UpdatedAtoms.dat"
"make_prt.dat")

for f in "${rm_files[@]}"; do
   rm -f -- $f 
done
