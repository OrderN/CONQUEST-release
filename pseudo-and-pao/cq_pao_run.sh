#!/bin/bash
#
# **<lat>** 
# Automated PAO generation for Conquest - 2020/07/15
#
# setup echo 
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
#
# setup env from file env_cq_pao
. ./src/env_cq_pao
#
CQ_ATOMLIST="cq_pao.list"
CQ_PAO_COMMAND=$CQ_BINDIR/MakeIonFiles
CQ_PAO_OUTPUT="Conquest_ion_out"
CQ_PAO_SRC=$CQ_PAODIR/src
CQ_ION_FILE="CQ.ion"
CQ_SPEC_FILE="CQ.spec"
CQ_BLOCK_FILE="CQ.block"
#
if [[ "$1" == "-h" || "$1" == ""  ]]; then
    $ECHO "  Usage: $ ./cq_pao_run.sh -xc dirname -l atom.list"
    $ECHO "          where dirname is the directory of the xc functional"
    $ECHO "                atom.list is the list of atoms"     
    $ECHO
    exit 1
elif [ "$1" == "-xc" ]; then
    
    if [ "$2" == "" ]; then
	$ECHO
	$ECHO "  ERROR: xc directory name is missing"
	$ECHO
	exit 1
    else
	XC_DIR=$2
	if [ ! -d "$XC_DIR" ]; then
	    $ECHO "  ERROR: xc directory" $XC_DIR "doesn't exists"
	    exit 1
	else
	    $ECHO "  xc directory:" $XC_DIR
	fi
    fi
else
    $ECHO
    $ECHO "  Check usage: $ ./cq_pao_run.sh -h"
    $ECHO
    exit 1
fi
#
if [[ "$3" == ""  ]]; then
    $ECHO "  Usage: $ ./cq_pao_run.sh -xc dirname -l atom.list"
    $ECHO "          where dirname is the directory of the xc functional"
    $ECHO "                atom.list is the list of atoms"     
    $ECHO
    exit 1

elif [ "$3" == "-l" ]; then
    
    if [ "$4" == "" ]; then
	$ECHO
	$ECHO "  WARNING: atom list is missing"
	$ECHO "  default:" $CQ_ATOMLIST "will be used"
	$ECHO	
    else
	CQ_ATOMLIST=$4
	$ECHO "  cq_atomlist:" $CQ_ATOMLIST
    fi
else
    $ECHO
    $ECHO "  Check usage: $ ./cq_pao_run.sh -h"
    $ECHO
    exit 1
fi
#
# Run MakeIonFiles for the set of atoms
#
while read file  ; do

    if [ -d "$XC_DIR/$file" ] ; then
	cd $XC_DIR/$file
	$ECHO
	$ECHO "  >> running the calculation for" $XC_DIR/$file "...\c"
	$CQ_PAO_COMMAND > $CQ_PAO_OUTPUT
	check_failure $?
	$ECHO "  done"
	$CQ_PAO_SRC/cq_generate_PAOblock.sh $file$CQ_ION_FILE
	cp $file$CQ_ION_FILE $CQ_PAODIR/lib
	cp $file$CQ_SPEC_FILE $CQ_PAODIR/lib
	cp $file$CQ_BLOCK_FILE $CQ_PAODIR/lib	
	cd ../..
    fi
	
done < $CQ_ATOMLIST
#
