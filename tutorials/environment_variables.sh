CQROOT='/Users/lioneltruflandier/CONQUEST-tutorials/bin'
export CQROOT
CQMPI='mpirun -np'
export CQMPI
CQNUM=1
export CQNUM
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
VMD_CMD='/Applications/VMD_2.0.0a7-pre2.app/Contents/MacOS/startup.command'
export VMD_CMD
VESTA_CMD='/Applications/VESTA/VESTA.app/Contents/MacOS/VESTA'
export VESTA_CMD

echo "number of process: $CQNUM"
echo "number of threads: $OMP_NUM_THREADS"

