#!/bin/bash
#
# **<lat>** 
# Automated G2 tests for Conquest - 2014/10/12
#
# setup echo 
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
#
# setup env from file env_cq_tests
. ./env_cq_g2
#
CQ_BENCHNAME=PBE
CQ_G2BASISKIND="CQ"
#
CQ_ROOT=`cd ../ ; pwd`
CQ_G2DIR=$CQ_ROOT/g2/g2_test
CQ_WORKDIR=$CQ_G2DIR/$CQ_BENCHNAME
CQ_POSTPDIR=$CQ_ROOT/g2
#
CQ_COMMAND="$MPI_PREFIX $MPI_NPROCS $CQ_BINDIR/Conquest"
#
CQ_BUILDG2FILE=build_g2_1997_neutral.sh
CQ_BUILDCQFILE=build_g2_cqinputs.sh
CQ_G2INPUTFILE=cq_g2input_param.sh
CQ_OUTFILELIST="InfoGlobal.dat\
	hilbert_make_blk.dat\ 
	supp_pao.dat\
        proc_atoms*\
	proc_partitions*\
        fort.*\
        *matrix.*\
        *chden*\
        input.log"
#
#
if [ ! -d $CQ_G2DIR ] ; then
   $ECHO "  $CQ_G2DIR doesn't exists: EXIT"
   exit 
else
   cd $CQ_G2DIR
fi
#
#
if [ -x $CQ_BUILDG2FILE ] ; then
   $ECHO
   $ECHO "  building the g2 input files"
   $CQ_G2DIR/$CQ_BUILDG2FILE
else
   $ECHO "  the build file $CQ_BUILDG2FILE is not executable or missing from $CQ_G2DIR: EXIT"
   exit
fi
#
#
$ECHO 
$ECHO "  *****   G2(1997) Test for Conquest   *****  "
$ECHO "                           v0.1 - 2014/10/20  "
$ECHO 
$ECHO "  Conquest benchname:    $CQ_BENCHNAME"
$ECHO "  working directory:     $CQ_WORKDIR"
$ECHO "  running Conquest as:   $CQ_COMMAND"
$ECHO  
$ECHO "  the G2(1997) list of neutral molecules is:"
#
# 
array=()
spec=()
while read file  ; do
#   echo $file
   array=($(<$file))
   nspe=${array[0]}
   let curs1=2+2*nspe
   let curs2=3+2*nspe
   let curs3=$curs2+1
   pola=${array[1]}
   spec=(${array[@]:2:$nspe})
   chag=${array[@]:curs1:1}
   spin=${array[@]:curs2:1}
   head=${array[@]:curs3}
   $ECHO "  $file: $head" 

done < g2_1997_neutral.list 
#
cd ..
#
if   [[ "$1" == "-l" && -n "$2" ]]; then
   g2list=$2
   $ECHO
   $ECHO "  you are running tests from the list: $g2list"
   $ECHO      

elif [[ "$1" == "-f" && -n "$2" ]]; then
   list=()
   list=("${@:2}")
   g2list=user.list
   for file in "${list[@]}" ; do
      $ECHO $file
   done > $g2list
   $ECHO
   $ECHO "  you are running tests from the following files:"
   $ECHO

elif [ "$1" == "-g2" ]; then
   g2list=g2_1997_neutral.list
   $ECHO
   $ECHO "  you are running tests from the full G2 list of neutral molecules"
   $ECHO

elif [ "$1" == "-h" ]; then
   $ECHO "  Usage:"
   $ECHO
   $ECHO "   $ check_cq_g2_neutral.sh -l file.list or -f file1 file2 file3..."
   $ECHO "     optional: you may provide"
   $ECHO "                               cq_g2input.param"
   $ECHO "                               basis_kind"
   $ECHO
   exit 1

else
   $ECHO
   $ECHO "  Check usage: $ check_cq_g2_neutral.sh -h"
   $ECHO
   exit 1

fi

if [ -f $CQ_G2INPUTFILE ]; then
   $ECHO "  user CQ input parameters file is:"
   $ECHO
   cat $CQ_G2INPUTFILE
   cp $CQ_G2INPUTFILE $CQ_G2DIR/cq_g2input_param.sh
   $ECHO
else
   $ECHO "  default CQ input parameters file will be used"
   $ECHO
   cp $CQ_G2DIR/cq_g2input_param_default.sh $CQ_G2DIR/cq_g2input_param.sh
   cat $CQ_G2DIR/cq_g2input_param.sh
   $ECHO
fi

#
#
if [ ! -f $g2list ] ;then
   cd $CQ_G2DIR 
   if [ ! -f $g2list ] ;then  
      $ECHO
      $ECHO "  ERROR: $g2list not existent"
      $ECHO "  Aborting"
      exit 1
   fi
else
   cp $g2list $CQ_G2DIR
   cd $CQ_G2DIR
fi

#
#
if [ ! -s $g2list ] ; then
   $ECHO
   $ECHO "  ERROR: $g2list is empty"
   $ECHO "  Aborting"
   exit 1
else
   while read line  ; do
      $ECHO "  $line"    
   done < $g2list
fi
#
#   
cd $CQ_G2DIR 
#
#
if [ -x $CQ_BUILDCQFILE ] ; then
   $ECHO
   $ECHO "  building the CQ input files"
   $CQ_G2DIR/build_g2_cqinputs.sh $g2list $CQ_G2INPUTFILE $CQ_G2BASISKIND
else
   $ECHO "  the build file $CQ_BUILDCQFILE is not executable or missing"
   $ECHO "  look at... $CQ_G2DIR"
fi

# With no arguments, checks all *.in files
# With an argument, checks files (ending with .in) matching the argument
#files=`/bin/ls *cq.in`
if [ ! "$CQ_BENCHNAME" == ""  ] ; then
   $ECHO
   $ECHO "  building benchmark directory"

   if [ ! -d "$CQ_BENCHNAME" ] ; then 
      mkdir $CQ_BENCHNAME
      $ECHO "  $CQ_WORKDIR has been created"
   else 
      $ECHO "  Warning: $CQ_WORKDIR working directory exists"
   fi
      while read file ; do 
         cqinput=`basename $file .dat`.cq.in
         cqcoord=`basename $file .dat`.coord

         mv $CQ_G2DIR/$cqcoord $CQ_WORKDIR
         mv $CQ_G2DIR/$cqinput $CQ_WORKDIR

      done < $g2list
 
      mv $CQ_G2DIR/*$CQ_G2BASISKIND.ion   $CQ_WORKDIR
      mv $CQ_G2DIR/*$CQ_G2BASISKIND.block $CQ_WORKDIR
      cp $CQ_G2DIR/$g2list $CQ_WORKDIR 

      $ECHO "  file(s) copied into $CQ_WORKDIR"
      cd $CQ_WORKDIR
      $ECHO "  now moved to $CQ_WORKDIR"
fi
#
$ECHO
#for file in $files
while read file
do
  name=`basename $file .dat`
  input=$name.cq.in
  output=$name.cq.out
  #
  if test -e Conquest_input ; then
      rm -f Conquest_input
  fi
  #
  ln -s $input Conquest_input
  #
  $ECHO "  Running $name..."
  ###
  # run the code in the G2 directory
  #
  #$MPI_PREFIX $MPI_NPROCS $CQ_ROOT/src/Conquest > $output 
  $CQ_COMMAND > $output
  #if test $? != 0; then
  #     $ECHO "  FAILED with error condition!"
  #     $ECHO "  Input: $name.in, Output: $name.out, Reference: $name.ref"
  #     $ECHO "  failed !"
  #     exit 1
  #else
  #   $ECHO " done"
  #fi
  #
  #cd $TESTDIR
  ###
  #if test -f $name.ref ; then
     # reference file exists
     #
     #check_scf $name
     #
     #get_times $name
     #get_nscf  $name
     #
  #else
  #   $ECHO  "not checked, reference file not available "
  #fi

  rm -f Conquest_input
done < $g2list

$ECHO
$ECHO "  Cleaning working directory: "
for file in $CQ_OUTFILELIST; do 
   if [ -f $file ]; then
      rm -f $file
      $ECHO "    $file removed"
   fi 
done


if [ ! "$CQ_BENCHNAME" == ""  ] ; then
   cd $CQ_G2DIR
   $ECHO
   $ECHO "  Archiving files for post-processing in: $CQ_BENCHNAME.tar  "
   tar -cf $CQ_BENCHNAME.tar $CQ_BENCHNAME
   mv $CQ_BENCHNAME.tar  $CQ_POSTPDIR
   $ECHO "  Archive moved to: $CQ_POSTPDIR "
fi

$ECHO
$ECHO "  Conquest benchname:    $CQ_BENCHNAME"
$ECHO "  working directory:     $CQ_WORKDIR"
$ECHO
$ECHO "  *****   G2(1997) Test for Conquest: $g2list done"
$ECHO

