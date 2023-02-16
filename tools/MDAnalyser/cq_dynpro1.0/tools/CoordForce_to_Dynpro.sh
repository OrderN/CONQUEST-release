#!/bin/bash
#

file=CoordForce.dat
nat=203

i=0
while read line ; do
    if (( $i == $nat )) ; then
       i=1
    else
       let i+=1
    fi

    if (( $i == 1 )) ; then
        echo $line | awk '{print $1, $2}' >> pos.dat
        echo $line | awk '{print $1, $2}' >> for.dat
        echo $line | awk '{print $1, $2}' >> vel.dat
        #echo $header
    fi
    #if (( $i == 2 )) ; then
    #    qspin=`echo $line | awk '{print $1, $2}'`
    #fi
    if (( $i > 1 )) ; then
    #    Z=`echo $line | awk '{print $1}'`
    #    xyz=`echo $line | awk '{printf "%.8f\n %.8f\n %.8f\n", $2, $3, $4}'`
    #    name=${chemsymb[$Z-1]}
    #    echo $name $xyz >> tmp
    #echo $name $xyz
    #fi
       echo $line | awk '{print $2, $3, $4  }' >> pos.dat
       echo $line | awk '{print $5, $6, $7  }' >> for.dat
       echo $line | awk '{print $8, $9, $10 }' >> vel.dat
       #echo $line
    fi 

 
done < $file

