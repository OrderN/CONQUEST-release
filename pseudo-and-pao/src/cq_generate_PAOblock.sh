#!/bin/bash

source ../../src/env_cq_pao 
source $CQ_PREFIX/pseudo-and-pao/src/cq_periodic_table.sh

function extrema() {
awk 'BEGIN {min=1000; max=0;}; \
{ if($1<min && $1 != "") min = $1;  \
  if($1>max && $1 != "") max = $1; }\
     END {print min, max}' $1
}

export extrema 

dr=0.2000
ds=0.2000
smax=8.0

for file in $1 ; do

   blockfile=`basename $file .ion`.block
   specifile=`basename $file .ion`.spec

   Z=`grep 'Atomic number'  $file | awk '{print $1}'`
   val=`grep 'Valence charge' $file | awk '{print $1}'`
   symb=${chemsymb[${Z%.*}-1]}
   mass=${chemmass[${Z%.*}-1]}

   Lnl=(); orb=(); cut=()
   Lnl=(`awk '/Lmax for basis/  {print $1, $2}' $file`)
   orb=(`awk '/orbital l, n, z/ {print $2, $1, $3, $4, $5}' $file`)
   cut=(`awk '/orbital l, n, z/ {x=NR+1;next}(NR<=x){print $3}' $file`)
   #Lnl=(`grep 'Lmax for basis'  $file | awk '{print $1, $2}'`)
   #orb=(`grep 'orbital l, n, z' $file | awk '{print $2, $1, $3, $4, $5}'`)

   #echo
   echo '  processing:' $blockfile
   echo '  Atom:' $symb '(Z='$Z'); lmax=' ${Lnl[0]}', nnl=' ${Lnl[1]}
   echo '  n l z p  occ   cutoff' 
   max=${Lnl[1]} ; j=0; norb=0 
   if test -f rcutoff ; then rm -f rcutoff ; fi
   for (( i = 0; i < max; i++ )); do 
      #echo 'j', $j, 'k', $k
      printf "  %d %d %d %d  %.2f  %.4f\n" ${orb[@]:$j:5} ${cut[$i]}
      echo ${cut[$i]} >> rcutoff 
      l=${orb[$j+1]}
      let "norb=$norb+(2*$l+1)"
      let j=$j+5
   done
   minmax=($(extrema rcutoff))
   min=${minmax[0]}
   max=${minmax[1]}
   max=`printf "%0.4f\n" $(bc -q <<< scale=4\;$max)`
   min=`printf "%0.4f\n" $(bc -q <<< scale=4\;$min)`
   echo "  norb=" $norb", rmin=" $min", rmax=" $max 
   rmax=`printf "%0.1f\n" $(bc -q <<< scale=4\;$max+$dr)`
   #smax=`printf "%0.4f\n" $(bc -q <<< scale=4\;$max+$ds)`

cat > $blockfile << EOF
%block $symb 
Atom.ValenceCharge    $val 
Atom.NumberOfSupports $norb 
Atom.SupportFunctionRange $rmax 
Atom.InvSRange            $smax 
%endblock 
EOF

cat > $specifile << EOF
%block ChemicalSpeciesLabel
 1  $mass  $symb  $file  
%endblock
EOF

done


