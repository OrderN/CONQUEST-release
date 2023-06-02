#!/bin/bash

# This script checks if the Conquest executable exists in ../bin and
# compiles it if it does not exist. Then it executes all tests found
# in directories named test_000_* to test_999_*, and runs pytest to
# check the results using test_check_output.py
#
# You can pass in the number of parallel processes as a command line
# argument, by default it is 1.

NP=${1:-1}
FILE=../bin/Conquest
if [ ! -f "$FILE" ]; then
    (cd ../src; make -j $NP)
fi

for dn in $(ls -d test_[0-9][0-9][0-9]_*)
do
    (cd $dn; mpirun -np $NP  ../../bin/Conquest)
done

pytest test_check_output.py
