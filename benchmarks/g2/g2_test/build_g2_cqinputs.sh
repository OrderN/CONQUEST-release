#!/bin/bash

if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
. ./../env_cq_g2

source convert_coord.sh
source periodic_table.sh
source cq_g2input_param.sh

if   [ "$1" == ""  ] ; then
   $ECHO
   $ECHO "  ERROR: g2list file is missing"
   $ECHO "  Aborting"
   exit 1
elif [ ! -s $1 ] ; then
   $ECHO
   $ECHO "  ERROR: g2list file is empty"
   $ECHO "  Aborting"
   exit 1
else
   g2list=$1
fi

if   [ "$2" == "" ] ; then
   $ECHO
   $ECHO "  ERROR: CQ input parameters file is missing"
   $ECHO "  Aborting"
   exit 1
elif [ ! -s $2 ] ; then
   $ECHO
   $ECHO "  ERROR: CQ input parameters file is empty"
   $ECHO "  Aborting"
   exit 1
else
   cqtemplate=$2
   cqchemspec=tmp_chemspec
fi

if   [ "$3" == "" ] ; then
   $ECHO
   $ECHO "  ERROR: CQ basis option is missing"
   $ECHO "  Aborting"
   exit 1
#elif [ ! -s $3 ] ; then
#   $ECHO
#   $ECHO "  ERROR: CQ basis file is empty"
#   $ECHO "  Aborting"
#   exit 1
else
   cqbasiskind=$3
fi
#
#
awk '{print $1}' $g2list > tmp 
mv tmp $g2list
#
#
array=()
spec=()
while read file  ; do
   # setup input and coordinate files
   cqinput=`basename $file .dat`.cq.in
   cqcoord=`basename $file .dat`.coord
   xyzfile=`basename $file .dat`.xyz
   # split data file to an array and extract
   #array=($(<$file))
   #nspe=${array[0]}
   #let curs1=2+2*nspe
   #let curs2=3+2*nspe
   #let curs3=$curs2+1
   #pola=${array[1]}
   #spec=(${array[@]:2:$nspe}) ; spec_sav=${spec[@]}
   #chag=${array[@]:curs1:1}
   #spin=${array[@]:curs2:1}
   #head=${array[@]:curs3}

   array=($(<$file))
   nspe=${array[0]}
   let curs1=2+2*nspe
   let curs2=3+2*nspe
   let curs3=4+2*nspe
   let curs4=$curs3+1
   pola=${array[1]}
   #nats=(${array[@]:2+$nspe:$nspe})
   spec=(${array[@]:2:$nspe}) ; spec_sav=${spec[@]}
   chag=${array[@]:curs1:1}
   spin=${array[@]:curs2:1}
   magn=${array[@]:curs3:1}
   head=${array[@]:curs4}

   header=$(printf ",%s" "${head[@]}")
   header=${header:1}

   # get nuclear charges z (needed to grab atom masses)
   z=() 
   for spe in ${spec[@]}; do
      i=0
      for symb in ${chemsymb[@]}; do
         if [ "$symb" ==  "$spe" ]; then
            z+=($i)         
         fi
      let i=$i+1
      done
   done
   # copy PAO ion/block files to working directory
   for spe in ${spec[@]}; do  
      basisfile=$spe$3.ion
      basisbloc=$spe$3.block 
      basispath=$CQ_PAODIR
      for paofile in "$basisfile" "$basisbloc"; do
         if [ ! -f "$basispath/$paofile" ] || [ ! -s "$basispath/$paofile" ]
         then
            $ECHO
            $ECHO "  ERROR: PAO  file(s) is missing or empty" 
            $ECHO "  Data file: $file"
            $ECHO "  Basis file: $basispath/$paofile"
            $ECHO "  Aborting"
            exit 1
            cp $basispath/$basisbloc . 
         else
            cp $basispath/$basisfile .
            cp $basispath/$basisbloc . 
         fi
      done 
   done
   # edit ChemicalSpeciesLabel block
   if test -f tmp_chemspec ; then rm -f tmp_chemspec ; fi
   echo                               >> $cqchemspec
   echo "%block ChemicalSpeciesLabel" >> $cqchemspec 
   i=1
   for spe in ${spec[@]}; do
      ztmp=${z[$i-1]}
      basisfile=$spe$3.ion
      #echo  "$i ${chemmass[$ztmp]} $spe  $basisfile"
      echo  "$i ${chemmass[$ztmp]} $spe  $basisfile" >> $cqchemspec 
      let i=$i+1  
   done 
   echo "%endblock" >> $cqchemspec 
   # grad atom block/ion files from basis directory
   

$ECHO "  CQ input: $cqinput =>\c" 
#$ECHO "  CQ input: $cqinput ..."
#
# edit Conquest input file
#cat > $cqinput << EOF
#IO.Title $head 
#IO.Coordinates $cqcoord  
#IO.FractionalAtomicCoords F
#IO.WriteOutToFile         F
#IO.Iprint                 1
#
## General Parameters
#General.NumberOfSpecies      $nspe 
#General.PseudopotentialType  siesta
#General.PartitionMethod      Hilbert
#General.PAOFromFiles         T
#
## Moving Atoms
#AtomMove.TypeOfRun static
#
## Spin Polarisation
#Spin.SpinPolarised F    
#Spin.FixSpin       F  
#Spin.NeUP          0.0 
#Spin.NeDN          0.0 
#
#EOF
#
#
# generate cq input file
#celldm=0.0
#cq_input celldm #$cqinput $cqcoord $nspe
#echo $header
if [ $spin ==  1 ]; then
   cq_input      celldm #$pwinput $nat $nspec $PW_TMPDIR $PW_PSEUDODIR
elif [ $spin -gt 1 ] ; then
   cq_input_spin celldm #$pwinput $nat $nspec $magn $PW_TMPDIR $PW_PSEUDODIR
else
   $ECHO
   $ECHO "  ERROR: CQ input spin parameters is wrong"
   $ECHO "  Aborting"
   exit 1
fi


#
# generate cq coord file 
xyz_to_coord $xyzfile $celldm
#
#cat $cqtemplate >> $cqinput  
cat $cqchemspec >> $cqinput

for spe in ${spec_sav[@]}; do
   basisbloc=$spe$3.block
#   echo $spe $basisbloc
   echo           >> $cqinput
   cat $basisbloc >> $cqinput 
done

#$ECHO "done"

done < $g2list 

if test -f tmp_chemspec ; then rm -f tmp_chemspec ; fi

