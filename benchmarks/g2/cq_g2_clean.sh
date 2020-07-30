#!/bin/sh

for dir in g2_test 
do
  cd $dir
  ./clean_g2_cqfiles.sh $1
  cd ..
done

