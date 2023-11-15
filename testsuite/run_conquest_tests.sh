#!/bin/bash

# This script checks if the Conquest executable exists in ../bin and
# compiles it if it does not exist. Then it executes all tests found
# in directories named test_000_* to test_999_*, and runs pytest to
# check the results using test_check_output.py
#
# You can pass in the number of parallel processes and number of
# OpenMP threads as a command line arguments, by default both are 1.

NP=${1:-1}
NT=${2:-1}
export OMP_NUM_THREADS=$NT
export OMP_STACKSIZE=100M
echo "Running tests on $NP processes and $NT threads"

(cd ../src; make -j $NP)

for dn in $(ls -d test_[0-9][0-9][0-9]_*)
do
    (cd $dn; mpirun -np $NP  ../../bin/Conquest)
done

pytest test_check_output.py
