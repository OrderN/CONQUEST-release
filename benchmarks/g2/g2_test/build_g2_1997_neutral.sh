#!/bin/bash

if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
. ./../env_cq_g2

source convert_coord.sh
source periodic_table.sh
source split_g2_1997_neutral.sh

g2filelist=g2_1997_neutral.list
#cbox=30.0
#	
# test if the 148 *.g2 files exist 
#ct=`ls -1 *.g2 2>/dev/null | wc -l`
ct=0
while read file ; do
file1=`basename $file .dat`.g2
  if [ -f $file1 ] ;then
    let ct+=1
  fi
done < $g2filelist
#
#
if [ ! $ct == 148 ] ; then
   $ECHO 
   $ECHO "  build_g2_1997_neutral/split_g2_neutral:" 
   $ECHO "  > files are missing or can be corrupted"
   $ECHO "  > $ct found / 148 expected"
   $ECHO 
   $ECHO "  splitting g2 raw data file..." 
   # split global g2 data file
   split_g2_neutral 
   $ECHO "  ...build_g2_1997_neutral/split_g2_neutral: done"
else
   $ECHO "  build_g2_1997_neutral/split_g2_neutral: nothing to be done"
fi
#
#
ct=0
while read file ; do
file1=`basename $file .dat`.xyz
  if [ -f $file1 ] ;then
    let ct+=1
  fi
done < $g2filelist
#
#
# test if the 148 *.xyz/ files exist 
#ct=`ls -1 *.xyz 2>/dev/null | wc -l`
if [ ! $ct == 148 ] ; then
   $ECHO
   $ECHO "  build_g2_1997_neutral/g2_to_xyz:"
   $ECHO "  > files are missing or can be corrupted"
   $ECHO "  > $ct found / 148 expected"
   $ECHO
   $ECHO "  converting g2 files to xyz..."

   while read file ; do
   #for g2file in *.g2 ; do
      xyzfile=`basename $file .dat`.xyz
      g2file=`basename $file .dat`.g2
      if [ ! -f $xyzfile ]; then
         g2_to_xyz $g2file
         #g2_to_cxyz $g2file $cbox # need to center the molecule within cell
      fi
   #done
   done < $g2filelist
   $ECHO "  ...build_g2_1997_neutral/g2_to_xyz: done"
else
   $ECHO "  build_g2_1997_neutral/g2_to_xyz: nothing to be done"
fi
#
#
ct=0
while read file ; do
  if [ -f $file ] ;then
    let ct+=1
  fi
done < $g2filelist
#
#
# test if the 296 *.coord + *.dat files exist 
#ct1=`ls -1 *.coord 2>/dev/null | wc -l`
#ct2=`ls -1 *.dat   2>/dev/null | wc -l`
#let "ct=$ct1+$ct2"
#let "ct=$ct2"
if [ ! $ct == 148 ] ; then
   $ECHO
   $ECHO "  build_g2_1997_neutral/xyz_to_dat: files are missing"
   $ECHO "  $ct found / 148 expected"
   $ECHO
   $ECHO "  converting xyz files to dat files..."

   while read file ; do
   #for g2file in *.g2 ; do
      #coordfile=`basename $g2file .dat`.coord
      datafile=`basename  $file .dat`.dat
      xyzfile=`basename   $file .dat`.xyz
      if [ ! -f $datafile ]; then
        xyz_to_dat $xyzfile
      fi
   #done
   done < $g2filelist
   $ECHO "  ...build_g2_1997_neutral/xyz_to_dat: done"
else
   $ECHO "  build_g2_1997_neutral/xyz_to_dat: nothing to be done"
fi

# define the list you want
# build a function to setup cq input
#   * need to read *.dat to get species / spin / charge
#   * define cards for species block
#   * ion files are within the working directory

#array=()
#spec=()
#while read file  ; do
#   echo $file
#   array=($(<$file))
#   nspe=${array[0]}
#   let curs1=2+2*nspe   
#   let curs2=3+2*nspe 
#
#   pola=${array[1]}
#   spec=(${array[@]:2:$nspe})
#   chag=${array[@]:curs1:1} 
#   spin=${array[@]:curs2:1}
#   echo $pola
#   echo ${spec[@]}
#   echo $chag
#   echo $spin

#done <g2_1997_neutral.list

check_failure $?
