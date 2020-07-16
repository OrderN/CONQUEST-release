#!/bin/bash
#
# **<lat>** 
# Automated PAO cleaning...
#
# setup echo #!/bin/bash
#
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
#
# setup env from file env_cq_pao
. ./src/env_cq_pao
#
CQ_ATOMLIST="cq_pao.list"
CQ_OUTFILELIST="Conquest_ion_out\
        input.log\
        rcutoff\
        *.block\
        *.spec\
        *.ion"
#
if [[ "$1" == "-h" || "$1" == ""  ]]; then
    $ECHO "  Usage: $ ./cq_pao_clean.sh -xc dirname"
    $ECHO "          where dirname is the directory of the xc functional"
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
    $ECHO "  Check usage: $ ./cq_pao_clean.sh -h"
    $ECHO
    exit 1
fi
$ECHO
$ECHO "  Cleaning PAO generation output data files:"
$ECHO
while read file  ; do
    #
    if [ -d "$XC_DIR/$file" ] ; then	
	cd $XC_DIR/$file
	for file1 in $CQ_OUTFILELIST; do
	    if [ -f $file1 ]; then
		rm -f $file1	   
		$ECHO "  $XC_DIR/$file/$file1 removed"
	    #else
		#pwd
		#$ECHO "  $XC_DIR/$file/$file1... nothing to do"
	    fi
	done
	cd ../..
    fi
    #
done < $CQ_ATOMLIST
#
cd lib
rm -f *.ion *.spec *.block
cd ..
