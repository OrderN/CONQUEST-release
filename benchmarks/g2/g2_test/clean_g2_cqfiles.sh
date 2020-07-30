#!/bin/bash
#
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
#
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


$ECHO
$ECHO "  Cleaning CQ output data files:"
$ECHO
for file in $CQ_OUTFILELIST; do
   if [ -f $file ]; then
      rm -f $file
      $ECHO "  $file removed"
   else 
      $ECHO "  $file... nothing to do"
   fi

done

if   [ "$1" == "medium" ]; then
   $ECHO
   $ECHO "  Cleaning CQ input/output/ion files:"
   $ECHO
   for file in *.cq.out *.cq.in *.ion *.block *.coord; do
      if [ -f $file ]; then
         rm -f $file
         $ECHO "  $file removed"
      else
         $ECHO "  $file... nothing to do"
      fi
   done

elif [ "$1" == "fine" ]; then
   $ECHO
   $ECHO "  Cleaning CQ input/output/ion/xyz/dat/coord files:"
   $ECHO
   for file in *.cq.out *.cq.in *.ion *.block *.coord *.cxyz *.xyz *.dat ; do
      if [ -f $file ]; then
         rm -f $file
         $ECHO "  $file removed"
      else
         $ECHO "  $file... nothing to do"
      fi
   done

elif [ "$1" == "full" ]; then
   $ECHO
   $ECHO "  Full cleaning:"
   $ECHO
   for file in *.cq.out *.cq.in *.ion *.block *.coord *.cxyz *.xyz *.dat *.g2; do
      if [ -f $file ]; then
         rm -f $file
         $ECHO "  $file removed"
      else
         $ECHO "  $file... nothing to do"
      fi
   done

elif [ "$1" == "light" ]; then
   $ECHO
   $ECHO "  Cleaning CQ input/output files:"
   $ECHO
   for file in *.cq.out *.cq.in ; do
      if [ -f $file ]; then
         rm -f $file
         $ECHO "  $file removed" 
      else
         $ECHO "  $file... nothing to do"
      fi
   done
else
   $ECHO "  nothing to do"
fi

$ECHO
$ECHO "  End cleaning"
$ECHO
