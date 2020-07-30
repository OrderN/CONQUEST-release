#!/bin/bash
#
# **< lat >** 16Oct2014
# 
source periodic_table.sh

#### function g2_to_xyz
#
# passed arguments:
#   *.g2  (string) g2 formated file
#
# return: 
#   *.xyz (string) xyz formated file
#
function g2_to_xyz() {
    xyzfile=`basename $1 .g2`.xyz
    i=0
    while read line ; do
	let i+=1
        if (( $i == 1 )) ; then
            header=`echo $line | awk '{$1=" "; print $0}'`
        fi
        if (( $i == 2 )) ; then
            qspin=`echo $line | awk '{print $1, $2}'`
        fi
	if (( $i > 2 )) ; then
	    Z=`echo $line | awk '{print $1}'`
	    xyz=`echo $line | awk '{printf "%.8f\n %.8f\n %.8f\n", $2, $3, $4}'`
	    name=${chemsymb[$Z-1]}
	    echo $name $xyz >> tmp
            #echo $name $xyz
	fi
    done < $1    
    let i=$i-2
    echo $i      >  $xyzfile
    echo $qspin $header >> $xyzfile
    cat tmp      >> $xyzfile
    tr -d '\b\r'  < $xyzfile > tmp
    mv tmp $xyzfile
    rm -rf tmp
}
#
#
#
#### function g2_to_cxyz
#
# passed arguments:
#   *.g2  (string) g2 formated file
#
# return: 
#   *.cxyz (string) cxyz formated file
#
function g2_to_cxyz() {
    cxyzfile=`basename $1 .g2`.cxyz

    size=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2)`
    cent=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2/2.0)`
    zero=`printf "%12.8f\n" $(bc -q <<< scale=8\;0.0)`

    i=0
    while read line ; do
        let i+=1
        if (( $i == 1 )) ; then
            header=`echo $line | awk '{$1=" "; print $0}'`
        fi
        if (( $i == 2 )) ; then
            qspin=`echo $line | awk '{print $1, $2}'`
        fi
        if (( $i > 2 )) ; then
            Z=`echo $line | awk '{print $1}'`
            #xyz=`echo $line | awk '{printf "%.8f\n %.8f\n %.8f\n", $2, $3, $4}'`
            name=${chemsymb[$Z-1]}

            x=`echo $line | awk '{printf "%.8f\n", $2}'`
            y=`echo $line | awk '{printf "%.8f\n", $3}'`
            z=`echo $line | awk '{printf "%.8f\n", $4}'`
            x=`printf "%0.8f\n" $(bc -q <<< scale=8\;$x+$cent)`
            y=`printf "%0.8f\n" $(bc -q <<< scale=8\;$y+$cent)`
            z=`printf "%0.8f\n" $(bc -q <<< scale=8\;$z+$cent)`

            echo $name $x $y $z >> tmp
            
            #echo $name $xyz >> tmp
            #echo $name $xyz
        fi
    done < $1
    #let i=$i-2
    #echo $i      >  $cxyzfile
    #echo $qspin $header >> $cxyzfile
    #cat tmp      >> $cxyzfile
    #tr -d '\b\r'  < $cxyzfile > tmp
    mv tmp $cxyzfile
    rm -rf tmp
}
#
#
#
#### function xyz_to_cq
#
# passed arguments:
#   cbox  (real kind) size of the cubic box  
#   *.xyz (string) name of the xyz file
#
# return:
#   *_cq.in  (string) Conquest input file
#   *_cq.xyz (string) Conquest coord file 
#
function xyz_to_coord() {
    bohr='0.52917726'
    cqcoord=`basename $1 .xyz`.coord
    #data=`basename $1 .xyz`.dat
    line=`head -1 $1` 
    nat=`echo $line | awk '{print $1}'`
    size=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2/$bohr)`
    cent=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2/$bohr/2.0)`
    zero=`printf "%12.8f\n" $(bc -q <<< scale=8\;0.0)`
    echo '   '$size'   ' $zero'   '$zero >  $cqcoord
    echo '   '$zero'   ' $size'   '$zero >> $cqcoord
    echo '   '$zero'   ' $zero'   '$size >> $cqcoord
    echo $nat >> $cqcoord 

    i=0 ; p=() ; qspin=()
    while read line ; do
        let i+=1
        #if (( $i == 2 )) ; then
        #    header=`echo $line | awk 'BEGIN {FS=" "; OFS=""} \
        #           {for(i=1;i<=NF;++i){out = out OFS $i}} END {print out;}'`    
        #    echo $header
        #fi 
        if (( $i == 2 )) ; then
            charge=`echo $line | awk '{print $1}'`
            spin=`echo $line   | awk '{print $2}'`
            header=`echo $line | awk '{print substr($0, index($0,$3))}'`
        fi 
        if (( $i > 2 )) ; then
            spec=`echo $line | awk '{print $1}'`
            p+=($spec)
        fi
    done < $1

    echo -e "  $header...\c"    

    q=() ; n=() ; s=() ; t=()
    #tail -$nat $2 | awk '{h[$1]++}; END { for(k in h) print k, h[k] }'
    n=(`tail -$nat $1 | awk '{h[$1]++}; END { for(k in h) print h[k] }'`)
    q=(`tail -$nat $1 | awk '{h[$1]++}; END { for(k in h) print k }'`)

    max=${#q[@]}
    for ((i = 1 ; i <= max ; i++ )); do s+=($i) ; done 
    #if (( $spin > 1)) ; then rad='U' ; else rad='R' ; fi
    #echo  $max $rad  > $data
    #echo  ${q[@]}   >> $data 
    #echo  ${n[@]}   >> $data
    #echo  $charge $spin >> $data
    #echo  $header       >> $data

    ct1=0  
    for i in ${q[@]} ; do
       ct2=0
       for j in ${p[@]} ; do
          if [ "$j" == "$i" ] ; then
              t[$ct2]=${s[$ct1]} 
          fi 
          let ct2+=1
       done
       let ct1+=1 
    done

    i=0
    while read line ; do
        let i+=1
        if (( $i > 2 )) ; then
            #xyz=`echo $line | awk '{printf "%.8f\n %.8f\n %.8f\n", $2, $3, $4}'`
            x=`echo $line | awk '{printf "%.8f\n", $2}'`
            y=`echo $line | awk '{printf "%.8f\n", $3}'`
            z=`echo $line | awk '{printf "%.8f\n", $4}'` 
            x=`printf "%0.8f\n" $(bc -q <<< scale=8\;$x/$bohr+$cent)`
            y=`printf "%0.8f\n" $(bc -q <<< scale=8\;$y/$bohr+$cent)`
            z=`printf "%0.8f\n" $(bc -q <<< scale=8\;$z/$bohr+$cent)` 
            name=${p[$i-3]}
            #echo $xyz ${t[$i-3]} 'T T T' >> $cqcoord 
            echo $x $y $z  ${t[$i-3]} ' T T T ' >> $cqcoord
            #echo "  " $name $xyz
        fi
    done < $1

    check_failure $? 
    echo " done"
}


function xyz_to_dat() {

    datafile=`basename $1 .xyz`.dat
    #datafile=$1

    bohr='0.52917726'
    line=`head -1 $1`
    nat=`echo $line | awk '{print $1}'`

    #size=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2/$bohr)`
    #cent=`printf "%12.8f\n" $(bc -q <<< scale=8\;$2/$bohr/2.0)`
    #zero=`printf "%12.8f\n" $(bc -q <<< scale=8\;0.0)`

    #echo '   '$size'   ' $zero'   '$zero >  $cqcoord
    #echo '   '$zero'   ' $size'   '$zero >> $cqcoord
    #echo '   '$zero'   ' $zero'   '$size >> $cqcoord
    #echo $nat >> $cqcoord

    i=0 ; p=() ; qspin=()
    while read line ; do
        let i+=1
        #if (( $i == 2 )) ; then
        #    header=`echo $line | awk 'BEGIN {FS=" "; OFS=""} \
        #           {for(i=1;i<=NF;++i){out = out OFS $i}} END {print out;}'`    
        #    echo $header
        #fi 
        if (( $i == 2 )) ; then
            charge=`echo $line | awk '{print $1}'`
            spin=`echo $line   | awk '{print $2}'`
            header=`echo $line | awk '{print substr($0, index($0,$3))}'`
        fi

        mag=`printf "%0.8f\n" $(bc -q <<< scale=8\;$spin-1.0)`

        if (( $i > 2 )) ; then
            spec=`echo $line | awk '{print $1}'`
            p+=($spec)
        fi
    done < $1

    echo -e "  $header...\c"    

    q=() ; n=() ; s=() ; t=()
    #tail -$nat $2 | awk '{h[$1]++}; END { for(k in h) print k, h[k] }'
    n=(`tail -$nat $1 | awk '{h[$1]++}; END { for(k in h) print h[k] }'`)
    q=(`tail -$nat $1 | awk '{h[$1]++}; END { for(k in h) print k }'`)

    max=${#q[@]}
    for ((i = 1 ; i <= max ; i++ )); do s+=($i) ; done
    if (( $spin > 1)) ; then rad='U' ; else rad='R' ; fi
    echo  $max $rad           > $datafile
    echo  ${q[@]}            >> $datafile
    echo  ${n[@]}            >> $datafile
    echo  $charge $spin $mag >> $datafile
    echo  $header            >> $datafile

    ct1=0
    for i in ${q[@]} ; do
       ct2=0
       for j in ${p[@]} ; do
          if [ "$j" == "$i" ] ; then
              t[$ct2]=${s[$ct1]}
          fi
          let ct2+=1
       done
       let ct1+=1
    done

    #i=0
    #while read line ; do
    #    let i+=1
    #    if (( $i > 2 )) ; then
    #       xyz=`echo $line | awk '{printf "%.8f\n %.8f\n %.8f\n", $2, $3, $4}'`
    #        x=`echo $line | awk '{printf "%.8f\n", $2}'`
    #        y=`echo $line | awk '{printf "%.8f\n", $3}'`
    #        z=`echo $line | awk '{printf "%.8f\n", $4}'`
    #        x=`printf "%0.8f\n" $(bc -q <<< scale=8\;$x/$bohr+$cent)`
    #        y=`printf "%0.8f\n" $(bc -q <<< scale=8\;$y/$bohr+$cent)`
    #        z=`printf "%0.8f\n" $(bc -q <<< scale=8\;$z/$bohr+$cent)`
    #        name=${p[$i-3]}
    #        #echo $xyz ${t[$i-3]} 'T T T' >> $cqcoord 
    #        echo $x $y $z  ${t[$i-3]} ' T T T ' >> $cqcoord
    #        #echo "  " $name $xyz
    #    fi
    #done < $1

    check_failure $?
    echo " done"
}






