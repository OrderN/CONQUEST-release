#!/bin/python
import os
import shutil
import subprocess
import time

subprocess.run(['bash', './../environment_variables.sh'])

# CQ root directory
CQROOT = os.getenv('CQROOT')
if (not CQROOT): 
    CQROOT='/Users/lioneltruflandier/CONQUEST-tutorials/bin'
# CQ mpi command
CQMPI = os.getenv('CQMPI')
if (not CQMPI): 
    CQMPI='mpirun -np'
# CQ nprocess
CQNUM = os.getenv('CQNUM')
if (not CQNUM): 
    CQNUM='2'
# CQ nthreads not greater than 4
OMP_NUM_THREADS = os.getenv('OMP_NUM_THREADS')
if ( not  OMP_NUM_THREADS ):
    os.environ['OMP_NUM_THREADS'] = '1'
# Setup Conquest command
CQCMD = CQMPI+' '+CQNUM+' '+CQROOT+'/Conquest'
# List of files to copy
files = ['Si.ion', 'coords.dat', 'Conquest_template']
# List of cutoff energy values in Ha
values=[50, 75, 100, 120, 150]
# basename
base = 'gs'
# Start timer
start = time.time()

print('Conquest cmd:',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname = str(val).rjust(3,'0')
    path = base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
        
    for file in files:
        shutil.copy2(file,path)

    os.chdir(path)    
    shutil.copy2('Conquest_template', 'Conquest_input')
    with open('Conquest_input','a') as f:
        f.write('\n\ngrid.gridcutoff '+str(val))
        
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')

# End timer
end = time.time()

# Print elapsed time
print('total of %.3f seconds elapsed'%(end -start))   
