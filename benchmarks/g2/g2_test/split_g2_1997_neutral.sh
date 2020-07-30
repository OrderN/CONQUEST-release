#!/bin/bash
#
function split_g2_neutral() {

linker='--link1--'
headerfile=header_g2test1997_neutral.rlist
structfile=struct_g2test1997_neutral.rlist
g2filelist=g2_1997_neutral.list

if [ ! -f $headerfile ] ; then
   echo 
   echo "ERROR: $headerfile file not existent or not in the directory"
   echo "Aborting"
   exit
fi

if [ ! -f $structfile ] ; then
   echo 
   echo "ERROR: $structfile file not existent or not in the directory"
   echo "Aborting"
   exit 1
fi

# split data G2/97 file
awk '/link1/{close(x);x=++i".g2";}{print > x;}' $structfile

# complete numerical basename
ct=0
while read file
do
  let ct+=1
  file1=$ct.g2
  file2=`echo "$file1" | awk -F_ '{ printf("%03d%s\n", $1, $2); }'`.g2
  #echo $file2
  if [ $file2 != $file1 ]; then
    mv $file1 $file2
  fi
done < $g2filelist

# add headers to structure files
ct=0
while read line ; do
  let ct+=1
  file=$(printf %03d.g2 ${ct%.g2})
  #echo "  $file" 
  sed -e "s/$linker/$linker $line/" $file > tmp ; mv tmp $file
done < $headerfile

}

