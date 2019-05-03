#!/bin/sh
#
# Adapted by Joao Lins, Instituto de Quimica da U.F.R.J.
#
# For really large computations, you can ultimately try:
# 1,\$s/numprm=400/numprm=1600/
#
# For more grid points try (For Gaussian and Vasp cube files)
# 1,\$s/mx3d=61/mx3d=122/
#
# To increase maximum allowed connections (In xwin.c change #define MXCON 6)
# 1,\$s/mxcon=6/mxcon=10/
#
# To increase maximum number of atoms for the Z-matrix (For converting
# proteins to Z-matrix) (Also change in xwin.c: #define MAXAT 1000)
# 1,\$s/maxat=1000/maxat=20000/
#
# 1,\$s/maxorb=256/maxorb=1024/
#
# IMPORTANT: molden3.6 has dynamic memory alloaction for the maxorb stuff
#            and for the protein Z-Matrix (the latter not for z-mat reading)
#            So the maxorb and maxat paramters should in general not be
#            tempered with
#
cat <<EOF > /tmp/sed.tmp
1,\$s/MAXPNT=1000/MAXPNT=200/
EOF
#
for file in *.f 
do
   echo Operating on $file
   mv $file ${file}.org
   sed -f /tmp/sed.tmp ${file}.org > $file
done
rm -f /tmp/sed.tmp
