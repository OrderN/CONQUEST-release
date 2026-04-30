#!/bin/bash

rm_files=("Conquest_out_ase" 
"Conquest_warnings"
"conquest.bib"
"coord_next.dat"
"eigenvalues.dat"
"hilbert_make_blk.dat"
"InfoGlobal.i00"
"input.log"
"Kmatrix*"
"UpdatedAtoms.dat"
"make_prt.dat"
"chden*"
"AtomCharge.dat"
"*cube"
"*WF*"
"Conquest_out_dos"
"Conquest_out_pp")

for f in "${rm_files[@]}"; do
   rm -f -- $f 
done
