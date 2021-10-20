#!/bin/bash
echo "Testing default PAO generation"
declare -i success=0
declare -i failure=0
declare -i total=0
for ele in Ag  As  Ba  Br  Cd  Cr  F   Ge  Hf  In  Kr  Mn  Na  Ni  P   Po  Re  Ru  Sc  Sn  Tc  Tl  Xe  Zr Al  Au  Be  C   Cl  Cs  Fe  H   Hg  Ir  Li  Mo  Nb  O   Pb  Pt  Rh  S   Se  Sr  Te  V   Y Ar  B   Bi  Ca  Co  Cu  Ga  He  I   K   Mg  N   Ne  Os  Pd  Rb  Rn  Sb  Si  Ta  Ti  W   Zn
do
  let total++
  cd $ele
  ../../../bin/MakeIonFiles > output 2>&1
  if [ `tail -1 output | grep -c Finish` -eq 1 ]; then
    let success++
    rm output "$ele"CQ.ion input.log
  else
    let failure++
    echo "Possible failure in default PAO generation for "$ele". Please report to developers."
    rm input.log
  fi
  cd ..
done
echo "Successes: "$success"/"$total
echo "Failures:  "$failure"/"$total
